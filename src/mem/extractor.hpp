#pragma once

#include "mem.hpp"
#include <string>

namespace lz77tax {

/**
 * Extractor de MEMs (Maximal Exact Matches) desde una query.
 *
 * Usa la BWT construida con ropebwt3 para extraer todos los MEMs
 * de una secuencia query contra el texto de referencia.
 */
class MEMExtractor {
public:
    /**
     * Extrae todos los MEMs de una query.
     *
     * @param query Secuencia a buscar
     * @param min_length Largo mínimo de MEM a reportar
     * @return Vector de MEMs encontrados
     */
    MEMList extract(const std::string& query, size_t min_length = 31);

private:
    // TODO: Agregar miembros para acceso a BWT
};

}  // namespace lz77tax
