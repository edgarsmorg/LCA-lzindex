/**
 * sr_locate — herramienta simple de locate sobre el sr-index.
 *
 * Uso:
 *   sr_locate --data_dir=<dir> --data_name=<basename> --patterns=<file>
 *
 * Formato de entrada (--patterns):
 *   Un patrón por línea (texto plano ASCII).
 *
 * Formato de salida (stdout):
 *   >PATRON
 *   pos1
 *   pos2
 *   ...
 *
 * Requiere que el índice ya esté construido en data_dir (archivos .sdsl
 * generados por bm_construct_ri).
 */

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include <gflags/gflags.h>
#include <sdsl/config.hpp>
#include <sdsl/int_vector_buffer.hpp>

#include "bench_common.hpp"
#include "sr-index/config.h"
#include "sr-index/r_index.h"

DEFINE_string(data_dir,  "./", "Directorio donde están los archivos .sdsl del índice.");
DEFINE_string(data_name, "",   "Nombre base del archivo de datos (OBLIGATORIO).");
DEFINE_string(patterns,  "",   "Archivo con patrones, uno por línea (OBLIGATORIO).");

int main(int argc, char** argv) {
    gflags::SetUsageMessage("Locate patterns in the sr-index.\n"
                            "Usage: sr_locate --data_dir=DIR --data_name=NAME --patterns=FILE");
    gflags::ParseCommandLineFlags(&argc, &argv, true);

    if (FLAGS_data_name.empty() || FLAGS_patterns.empty()) {
        std::cerr << "Error: --data_name y --patterns son obligatorios.\n";
        return 1;
    }

    // Cargar R-Index desde archivos de caché sdsl
    sri::Config config(FLAGS_data_name, FLAGS_data_dir,
                       sri::SDSL_LIBDIVSUFSORT, /*load=*/true, /*width=*/8);

    using ExternalStorage = std::reference_wrapper<sri::GenericStorage>;
    sri::GenericStorage storage;
    sri::RIndex<ExternalStorage> ri(std::ref(storage));

    try {
        ri.load(config);
    } catch (const std::exception& e) {
        std::cerr << "Error cargando el índice: " << e.what() << "\n";
        std::cerr << "Asegúrate de haber corrido bm_construct_ri primero.\n";
        return 2;
    }

    // Leer patrones y hacer locate
    std::vector<std::string> patterns;
    try { patterns = bench::read_patterns(FLAGS_patterns); }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 3;
    }

    using Sequence = sri::Alphabet<>::string_type;

    for (const auto& line : patterns) {
        // Convertir string a Sequence (uint8_t)
        Sequence pattern(line.begin(), line.end());

        auto positions = ri.Locate(pattern);
        std::sort(positions.begin(), positions.end());

        std::cout << ">" << line << "\n";
        for (auto pos : positions) {
            std::cout << pos << "\n";
        }
    }

    return 0;
}
