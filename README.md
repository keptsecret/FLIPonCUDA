# FLIP on CUDA

This project implements PIC/FLIP fluid simulation on CUDA/C++17.
The simulation program is intended to simulate liquids and outputs the surface of the fluid
as a sequence of triangles meshes.

The program features GPU acceleration on Nvidia cards supporting CUDA.
The two GPU-accelerated simulation steps are velocity field advection
and marker particle advection using RK4 integration.

## Build instructions

The dependencies required for the project are:

* C++17 compiler
* Nvidia CUDA Toolkit
  * this project was developed on version 11.7

The project uses CMake to generate the project files.
By default, the project is configured to build with GPU acceleration.
You can decide to build without GPU acceleration and run single-threaded on CPU by setting the CMake option:

```bash
cmake -DCUDA_ENABLED:BOOL=OFF
```

## Usage

You can run the program through the generated executable.
By default, the program saves the meshes to a `cache` folder inside the build directory.
Make sure to create the directory before running the executable.

```bash
mkdir cache
```

## Samples

![water sphere drop](./samples/flipsphere.gif)
*Simulated on 64x64x64 grid with cell size 0.125 for 60 frames at 30 FPS. Rendered in Blender Cycles.*
