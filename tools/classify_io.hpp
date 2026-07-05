#pragma once

/**
 * classify_io.hpp — loaders TSV compartidos por las herramientas de clasificación
 * (compare_classifiers, bench_classify_compare).
 *
 * Formatos (generados por scripts/tree_json_to_tsv.py o gen_bench_datasets.py):
 *   tree.tsv     node_id  parent_id  name  genome_idx
 *   genomes.tsv  leaf_dfs_rank  start  end  node_id
 *   reads.tsv    read_id  true_genome_idx  true_node_id  sequence
 */

#include "taxonomy/lca.hpp"
#include "taxonomy/classifier.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace classify_io {

/// Carga el árbol filogenético desde tree.tsv. Devuelve false si no puede abrirlo.
inline bool load_tree(const std::string& path, lz77tax::PhyloTree& tree) {
    std::ifstream f(path);
    if (!f) { std::cerr << "No se puede abrir: " << path << "\n"; return false; }

    std::vector<int>         ids, parents;
    std::vector<std::string> names;

    std::string line;
    std::getline(f, line);  // header
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string nid_s, pid_s, name, gidx_s;
        std::getline(ss, nid_s, '\t');
        std::getline(ss, pid_s, '\t');
        std::getline(ss, name,   '\t');
        std::getline(ss, gidx_s, '\t');
        ids.push_back(std::stoi(nid_s));
        parents.push_back(std::stoi(pid_s));
        names.push_back(name);
    }
    tree.build(ids, names, parents);
    return true;
}

/// Carga los rangos de genoma desde genomes.tsv.
inline bool load_genomes(const std::string& path,
                         std::vector<lz77tax::Classifier::GenomeRange>& genomes) {
    std::ifstream f(path);
    if (!f) { std::cerr << "No se puede abrir: " << path << "\n"; return false; }

    std::string line;
    std::getline(f, line);  // header
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string rank_s, start_s, end_s, nid_s;
        std::getline(ss, rank_s,  '\t');
        std::getline(ss, start_s, '\t');
        std::getline(ss, end_s,   '\t');
        std::getline(ss, nid_s,   '\t');
        lz77tax::Classifier::GenomeRange g;
        g.start   = std::stoull(start_s);
        g.end     = std::stoull(end_s);
        g.node_id = std::stoi(nid_s);
        genomes.push_back(g);
    }
    return true;
}

/// Una read con ground-truth taxonómico.
struct ReadEntry {
    int         read_id;
    int         true_genome_idx;
    int         true_node_id;
    std::string sequence;
};

/// Carga reads.tsv con ground-truth.
inline bool load_reads(const std::string& path, std::vector<ReadEntry>& reads) {
    std::ifstream f(path);
    if (!f) { std::cerr << "No se puede abrir: " << path << "\n"; return false; }

    std::string line;
    std::getline(f, line);  // header
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string id_s, gidx_s, nid_s, seq;
        std::getline(ss, id_s,   '\t');
        std::getline(ss, gidx_s, '\t');
        std::getline(ss, nid_s,  '\t');
        std::getline(ss, seq,    '\t');
        ReadEntry e;
        e.read_id         = std::stoi(id_s);
        e.true_genome_idx = std::stoi(gidx_s);
        e.true_node_id    = std::stoi(nid_s);
        e.sequence        = seq;
        reads.push_back(std::move(e));
    }
    return true;
}

}  // namespace classify_io
