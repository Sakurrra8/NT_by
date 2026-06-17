# Local Testing

This repository is still primarily a Linux/HPC codebase. The lowest-friction local setup is WSL2 Ubuntu, not native Windows.

## Python Triangle Benchmark

From the repository root:

```bash
python3 triangle/EAST/benchmark_tri.py --output Outputfile/test
```

The script now defaults `--root` to the current directory and writes figures to:

```text
Outputfile/test/fig/tri_benchmark/
```

You need output data copied from the HPC run, for example:

```text
Outputfile/test/data/
```

## C++ Build Under WSL/Linux

The normal `Makefile` is kept for the HPC environment. For local Linux or WSL builds, use:

```bash
make -f Makefile.local -j
```

Debug build:

```bash
make -f Makefile.local BUILD=debug -j
```

The repository includes a Linux-style HDF5 tree under `opt/hdf5`, so the default is:

```text
HDF5_DIR=./opt/hdf5
```

On Windows clones, Git may materialize Linux symlinks such as `libhdf5.so` as tiny text files. `Makefile.local` handles this by creating WSL-only symlinks under:

```text
.local/hdf5/lib/
```

If you use system HDF5 instead:

```bash
sudo apt install build-essential make libgomp1 libhdf5-dev
make -f Makefile.local HDF5_DIR=/usr -j
```

The executable is written to:

```text
bin/edit
```

Run the default triangle test case:

```bash
make -f Makefile.local run-check
```

Use another setting file:

```bash
make -f Makefile.local run-check SETTING=Inputfile/settingfile/setting_D.log
```

## VS Code Debugging With WSL

Install these VS Code extensions:

```text
WSL
C/C++
```

In Ubuntu/WSL:

```bash
sudo apt update
sudo apt install build-essential gdb make libgomp1
```

Open the repository through the VS Code WSL window:

```text
/mnt/d/Users/b11/Documents/codex/nt
```

The repository includes:

```text
.vscode/tasks.json
.vscode/launch.json
```

Press `F5` and choose:

```text
Debug nt
```

This builds with:

```bash
make -f Makefile.local BUILD=debug -j
```

and launches:

```bash
./bin/edit Inputfile/settingfile/setting_Trimesh_D_5.log
```

## Native Windows

Native Windows compilation is not the recommended path because the code includes Linux headers such as `unistd.h`, uses OpenMP flags from `g++`, and links Linux-style HDF5 libraries. Use WSL for local execution.
