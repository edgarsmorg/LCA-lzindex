# Memoria: Mapeo de secuencias de ADN a árboles filogenéticos usando compresión Lempel-Ziv

Proyecto de Memoria de Edgar S. Morales González para optar al título de Ingeniero Civil en Computación, Universidad de Chile, DCC.

**Profesor Guía:** Gonzalo Navarro

## Descripción

Construcción de un **índice comprimido basado en LZ77** para clasificación taxonómica de lecturas de ADN. El índice utiliza:

- Parsing LZ77 del texto concatenado
- Grilla bidimensional de puntos primarios
- Wavelet Tree con RMQ sucinto integrado
- Búsqueda de MEMs y mapeo a árboles filogenéticos

## Compilación

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel $(nproc)
```

## Dependencias

- **C++17** o superior
- **GCC 12+**
- sdsl-lite (fork duscob)
- ropebwt3
- sr-index (para comparación)

Véase [CLAUDE.md](CLAUDE.md) para detalles completos de configuración y arquitectura.
