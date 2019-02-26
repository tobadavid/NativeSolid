# NativeSolid
Native solid solver for FSI

## Features
Compute static or dynamic displacements for several 2D configurations.

* Meshes
  - .su2 format (can be built and exported from gmsh)
* Configuration
  - Airfoil (pitch/plunge)
  - Horizontal spring
  - Vertical spring
* Dynamic computation
  - AlphaGen time integration
  - Runge-Kutta (order 4) time integration
  - Mass/damping/stiffness/nonlinear terms
* Static Computation
  - Stiffness term only

## Compilation (linux -gcc)
Installation of NativeSolid code

Required packages
```
sudo apt-get install build-essential
sudo apt-get install cmake
sudo apt-get install liblapacke-dev
sudo apt-get install liblatlas-base-dev
```
Optional packages (parallel build)
```
sudo apt-get install openmpi-bin openmpi-common libopenmpi-dev
```
Compilation (you might need to set paths for headers and libraries in the FindCBLAS and FindLAPACKE files)
```
mkdir build && cd build
cmake ..
make -j4
```
Environment variables
```
gedit /home/"username"/.bashrc
export NATIVE_RUN="INSTALL_PATH/NativeSolid/bin"
export PATH=$PATH:$NATIVE_RUN
export PYTHONPATH=$PYTHONPATH:$NATIVE_RUN
	
