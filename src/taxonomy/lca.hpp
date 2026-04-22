#pragma once

#include <vector>
#include <string>
#include <cassert>

namespace lz77tax {

/**
 * Árbol filogenético con cómputo de LCA mediante subida naive.
 *
 * Para el prototipo sintético se usa subida por niveles (O(depth)).
 * En producción se reemplaza por LCA de Farach-Colton & Bender (O(1)).
 */
class PhyloTree {
public:
    static constexpr int NO_PARENT = -1;

    struct Node {
        int         id;
        std::string name;
        int         parent;
        int         depth;
    };

    PhyloTree() = default;

    /**
     * Construye el árbol desde listas paralelas de ids, nombres y padres.
     * El nodo raíz debe tener parent = NO_PARENT.
     */
    // REQUISITO: ids debe ser {0, 1, 2, ..., n-1} en el mismo orden que el vector.
    // lca() usa node_id directamente como índice de array; si ids[i] != i el
    // comportamiento es indefinido. parents[i] también debe ser un id válido (o NO_PARENT).
    void build(const std::vector<int>& ids,
               const std::vector<std::string>& names,
               const std::vector<int>& parents) {
        assert(ids.size() == names.size() && ids.size() == parents.size());
        int n = ids.size();
        nodes_.resize(n);
        for (int i = 0; i < n; ++i) {
            assert(ids[i] == i && "PhyloTree: ids deben ser secuenciales 0..n-1");
            nodes_[i] = {ids[i], names[i], parents[i], -1};
        }
        // Calcular profundidades desde la raíz
        for (auto& nd : nodes_) {
            if (nd.parent == NO_PARENT) {
                nd.depth = 0;
            }
        }
        bool changed = true;
        while (changed) {
            changed = false;
            for (auto& nd : nodes_) {
                if (nd.depth == -1 && nodes_[nd.parent].depth != -1) {
                    nd.depth = nodes_[nd.parent].depth + 1;
                    changed = true;
                }
            }
        }
    }

    /** LCA de dos nodos (subida naive por nivel). */
    int lca(int u, int v) const {
        // Subir al mismo nivel
        while (nodes_[u].depth > nodes_[v].depth) u = nodes_[u].parent;
        while (nodes_[v].depth > nodes_[u].depth) v = nodes_[v].parent;
        // Subir juntos hasta encontrar el ancestro común
        while (u != v) {
            u = nodes_[u].parent;
            v = nodes_[v].parent;
        }
        return u;
    }

    const std::string& name(int node_id) const { return nodes_[node_id].name; }
    int depth(int node_id) const { return nodes_[node_id].depth; }
    int size() const { return nodes_.size(); }

private:
    std::vector<Node> nodes_;
};

}  // namespace lz77tax
