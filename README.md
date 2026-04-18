# BoneSimulation (C++)

A C++20 processing pipeline for micro-CT image stacks of bone specimens.
The project loads TIFF stacks, binarises the volume (Otsu), performs
morphological operations and connected-component analysis, generates
surface and tetrahedral meshes, and exports the results for subsequent
finite element analysis with **CalculiX**.

The entire build and execution workflow runs inside a **Docker container**,
so no dependencies need to be installed on the host system.

---

## Table of Contents

- [Features](#features)
- [Quick Start](#quick-start)
- [Test Modes](#test-modes)
- [Project Layout](#project-layout)
- [Pipeline Overview](#pipeline-overview)
- [Dependencies](#dependencies)
- [Docker](#docker)
- [Local Build (optional)](#local-build-optional)
- [Input and Output Directories](#input-and-output-directories)
- [Command-Line Options](#command-line-options)
- [Development](#development)
- [Troubleshooting](#troubleshooting)
- [License](#license)

---

## Features

- **TIFF stack loader** (16-bit) for both directories of slices and multi-page TIFF files
- **Otsu thresholding** (sample-based) with in-place binarisation
- **Morphological closing** (3D reference implementation and slice-optimised variant)
- **Volume interpolation** (direct and chunked, for large volumes)
- **Connected-component analysis** (Union-Find, 26-neighbourhood)
- **Marching Cubes** surface extraction (VTK) with optional smoothing
- **Tetrahedral meshing** via `tetgen`
- **Parallel block extraction** (OpenMP) for sub-volume processing
- **CalculiX export** (`.inp` input deck plus runner script)
- Timestamped logging utility
- Fully containerised workflow (Docker and docker-compose)

---

## Quick Start

```bash
# 1) Build the image
docker build -f docker/Dockerfile -t bonesim:latest .

# 2) Run the full pipeline
docker run --rm \
    -v $(pwd)/data:/workspace/data \
    -v $(pwd)/pictures:/workspace/pictures \
    bonesim:latest --test complete
```

Input TIFF data is expected under `./data/bigdataset/`. Generated artifacts
(meshes, extracted blocks, etc.) will appear in `./pictures/`.

A convenience wrapper is also provided:

```bash
./scripts/run_docker.sh --test complete
```

---

## Test Modes

The `bonesim` binary accepts several test modes through the `--test <mode>` flag:

| Mode       | Description                                                                  |
|------------|------------------------------------------------------------------------------|
| `dry`      | Dry run; no files are written                                                |
| `small`    | Reduced pipeline (interpolation factor 0.2, significantly faster)            |
| `complete` | Full pipeline including CalculiX export                                      |

Example — quick smoke test:

```bash
docker run --rm \
    -v $(pwd)/data:/workspace/data \
    -v $(pwd)/pictures:/workspace/pictures \
    bonesim:latest --test small
```

---

## Project Layout

```text
BoneSimulation_cpp/
├── CMakeLists.txt              # Top-level CMake (delegates to src/)
├── README.md
├── docker/
│   ├── Dockerfile              # Production image (Ubuntu 22.04, multi-stage)
│   ├── Dockerfile.dev          # Development image (Arch Linux, full toolchain)
│   └── docker-compose.yml      # Compose setup including all volumes
├── scripts/
│   ├── build.sh                # Local release build (cmake + make)
│   ├── run_docker.sh           # Build and run in one step
│   └── install_deps_arch.sh    # Host dependency installer (Arch Linux)
├── src/
│   ├── CMakeLists.txt          # Actual `bonesim` target definition
│   ├── main.cpp                # Entry point, argument parsing, pipeline driver
│   ├── core/                   # Core algorithms
│   │   ├── image_loading.*     # TIFF stack / volume loaders
│   │   ├── image_processing.*  # Otsu, closing, interpolation, cropping
│   │   ├── cluster_analysis.*  # Union-Find connected components
│   │   ├── mesh_generation.*   # Marching Cubes, tetgen wrapper, VTK I/O
│   │   ├── block_extraction.*  # Parallel sub-volume extraction
│   │   ├── calculix_export.*   # `.inp` export and runner script generation
│   │   └── reporting.*         # Optional PDF report (libharu)
│   └── utils/
│       ├── file_io.*           # TIFF writer
│       └── logger.*            # Lightweight logger
├── data/                       # Input data (bind-mounted into the container)
│   └── bigdataset/             # TIFF stack (*.tif / *.tiff)
├── pictures/                   # Output: meshes, extracted blocks, intermediates
├── report/                     # Output: optional PDF reports
├── calculix/                   # Output: `.inp` decks and runner scripts
└── tests/                      # Reserved for future unit/integration tests
```

---

## Pipeline Overview

```text
TIFF stack (uint16)
    │
    ▼
Otsu threshold  ──►  in-place binarisation (→ uint8)
    │
    ▼
Morphological closing (slice-based)
    │
    ▼
Interpolation (chunked, configurable factor)
    │
    ▼
Largest connected component (Union-Find)
    │
    ├──► Marching Cubes  ───────────► pictures/mesh_<ts>.vtk
    │
    ├──► Block extraction (OpenMP) ─► pictures/blocks_<ts>/block_*.npy
    │                                   + block_*_tetra_<ts>.vtk
    │
    └──► Tetrahedral mesh ─────────►  calculix/bone_<ts>.inp
                                      calculix/run_simulation.sh
```

---

## Dependencies

All dependencies are provided by the Docker image:

- **Compiler:** `g++` (C++17/20)
- **Build system:** CMake ≥ 3.20, Make
- **Libraries:**
    - [VTK 9](https://vtk.org/) — `CommonCore`, `CommonDataModel`, `FiltersCore`, `IOGeometry`, `IOLegacy`
    - [Eigen 3](https://eigen.tuxfamily.org/)
    - [libtiff](http://www.libtiff.org/)
    - [Boost](https://www.boost.org/) — `system`, `thread`
    - [OpenMP](https://www.openmp.org/)
    - [tetgen](https://wias-berlin.de/software/tetgen/) (invoked as an external CLI)
    - *(optional)* [libharu](http://libharu.org/) for PDF report generation

---

## Docker

### Production image (`docker/Dockerfile`)

A multi-stage build based on **Ubuntu 22.04**. The build stage compiles the
binary; the runtime stage contains only the required shared libraries.

```bash
docker build -f docker/Dockerfile -t bonesim:latest .
```

Run:

```bash
docker run --rm \
    -v $(pwd)/data:/workspace/data \
    -v $(pwd)/pictures:/workspace/pictures \
    bonesim:latest --test complete
```

### docker-compose

A preconfigured setup with all volumes is available under `docker/docker-compose.yml`:

```bash
cd docker
docker compose up --build
```

### Development image (`docker/Dockerfile.dev`)

An Arch-based image including `gdb`, `valgrind`, `clang`, etc. for
interactive development:

```bash
docker build -f docker/Dockerfile.dev -t bonesim:dev .
docker run --rm -it -v $(pwd):/workspace bonesim:dev
```

---

## Local Build (optional)

If you prefer to build outside of Docker:

```bash
# Install dependencies (Arch Linux)
./scripts/install_deps_arch.sh

# Release build
./scripts/build.sh
# → Binary will be placed at ./build/bonesim
```

On Ubuntu/Debian the equivalent packages are those listed in
`docker/Dockerfile` (`build-essential cmake libtiff-dev libeigen3-dev
libvtk9-dev libboost-all-dev libopenmpi-dev tetgen`).

---

## Input and Output Directories

Inside the container the working directory is `/workspace`. The following
host paths are bind-mounted via `-v`:

| Host path    | Container path         | Purpose                             |
|--------------|------------------------|-------------------------------------|
| `./data`     | `/workspace/data`      | Input TIFFs (`data/bigdataset`)     |
| `./pictures` | `/workspace/pictures`  | Meshes and extracted blocks         |
| `./report`   | `/workspace/report`    | Optional PDF reports                |
| `./calculix` | `/workspace/calculix`  | CalculiX `.inp` decks and scripts   |

Generated files are tagged with a Unix timestamp, e.g.
`pictures/mesh_1713441200.vtk`.

---

## Command-Line Options

```text
Usage: bonesim [options]
  --test dry             Dry run (no files are written)
  --test small           Small test run (interpolation factor 0.2)
  --test complete        Full pipeline including CalculiX export
  --no-blocks            Skip block extraction
  --no-calculix          Skip CalculiX export
  --threads N            Set the number of OpenMP threads
  --help                 Show help
```

Default paths (relative to the container working directory):

- `data/bigdataset` — input
- `pictures`        — mesh output
- `report`          — PDF reports
- `calculix`        — FEM input

---

## Development

The project is a standard CMake project and can be opened in **CLion**,
**VS Code**, or any IDE with CMake support. The language standard is set
to **C++17** (`CMAKE_CXX_STANDARD 17`); the surrounding solution uses
C++20 features.

Typical edit/debug cycle:

```bash
# 1) Start a development container
docker build -f docker/Dockerfile.dev -t bonesim:dev .
docker run --rm -it -v $(pwd):/workspace bonesim:dev

# 2) Build inside the container
./scripts/build.sh

# 3) Run
./build/bonesim --test small
```

The `tests/` directory is reserved for future automated tests.

---

## Troubleshooting

**`No TIFF files in data/bigdataset`**
Ensure that `./data/bigdataset/*.tif` exists on the host and that the
directory is correctly bind-mounted into the container.

**`tetgen failed with code …`**
The surface mesh is most likely too large or not watertight. Retry with
`--test small` or a smaller interpolation factor.

**High memory usage**
The pipeline aggressively releases intermediate volumes, but very large
datasets may still exceed available RAM. In that case use `--test small`,
or disable stages via `--no-blocks` / `--no-calculix`.

**OpenMP uses only one core**
Set `--threads N` explicitly or pass `OMP_NUM_THREADS=N` as an
environment variable.

---

## License

*To be defined.* Please add an appropriate license file (e.g. `LICENSE`)
before publishing or distributing the project.
