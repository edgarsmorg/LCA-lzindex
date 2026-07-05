#pragma once

#include <climits>
#include <cstddef>
#include <memory>
#include <string>
#include <utility>

namespace lz77tax {

/**
 * Wrapper del sr-index (SrIndexValidArea) que expone la misma API de
 * locate_extremal que LZ77Index.
 *
 * Estrategia: ri.Locate(P) enumera todas las ocurrencias de P via Phi backward;
 * el mínimo y máximo se obtienen como min/max sobre el vector resultante.
 * Complejidad: O(m · r/s) por query donde r = runs de la BWT, s = subsampling rate.
 *
 * Usa PIMPL para aislar los headers del sr-index en sr_index_locator.cpp,
 * evitando colisiones de definición con index.cpp cuando se linkean juntos.
 */
class SrIndexLocator {
 public:
    SrIndexLocator();
    ~SrIndexLocator();

    SrIndexLocator(SrIndexLocator&&) noexcept;
    SrIndexLocator& operator=(SrIndexLocator&&) noexcept;

    SrIndexLocator(const SrIndexLocator&)            = delete;
    SrIndexLocator& operator=(const SrIndexLocator&) = delete;

    /**
     * Carga el índice desde archivos pre-construidos en data_dir.
     * Los archivos deben haber sido generados por sri::construct (o bm_construct_ri).
     *
     * @param data_name  Nombre base del dataset (e.g. "scaling_10MB")
     * @param data_dir   Directorio con los archivos .sdsl
     * @param sr         Subsampling rate usado al construir (default: 16)
     */
    void load(const std::string& data_name,
              const std::string& data_dir,
              std::size_t sr = 16);

    /**
     * Construye el índice desde un archivo de texto plano y lo mantiene en memoria.
     * Útil para tests sobre textos pequeños/sintéticos.
     *
     * @param data_path  Ruta al archivo de texto (texto plano, no FASTA)
     * @param cache_dir  Directorio temporal para archivos sdsl intermedios
     * @param sr         Subsampling rate (default: 16)
     */
    void build(const std::string& data_path,
               const std::string& cache_dir,
               std::size_t sr = 16);

    /**
     * Retorna {pos_min, pos_max} de todas las ocurrencias del patrón.
     *
     * Idéntica semántica a LZ77Index::locate_extremal, pero sin restricción
     * primary-only: ve todas las ocurrencias (primarias + secundarias).
     * Retorna {SIZE_MAX, 0} si el patrón no aparece en el texto.
     */
    std::pair<std::size_t, std::size_t>
    locate_extremal(const std::string& pattern) const;

    /**
     * Número de ocurrencias del patrón.
     * Útil para comparar con LZ77Index::count() (que cuenta solo primarias).
     */
    std::size_t locate_count(const std::string& pattern) const;

    /** Tamaño serializado del índice en bytes (footprint en disco). */
    std::size_t size_in_bytes() const;

    bool is_loaded() const;

 private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace lz77tax
