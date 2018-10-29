# PeriVIS

## Paraview plugins

### fragment 

Extraction of fragments and histograms, see Section 4.1 in [1]. The Figure 1 shows the extracted fracments and Figure 2 shows the corresponding histogram. A video of the development of the extracted fragments is available [here](https://archive.ins.uni-bonn.de/schweitzer.ins.uni-bonn.de/people/diehl/video/fragments.mp4).

Figure 1: Extracted fragments | Figure 2: Corresponding histogram
:------------------------:|:----------------------------------------:
![Fragments](./figures/fragments_20.jpg?raw=true "Fragments")|![Histogram](./figures/histo_19.jpg?raw=true "Histogram")


```
-DBUILD_FRAGMENT=ON 
```

### stress 

Volume rendering of the nodes and visualizaiton of the vonMises stress tensor, see Section 4.2 in [1]. Figures shows the volume rendering of the Standford bunny and the vonMises stress. An example video is available [here](https://archive.ins.uni-bonn.de/schweitzer.ins.uni-bonn.de/people/diehl/video/bunny.mp4).


Figure 1: Extracted fragments | 
:------------------------:|
![Stress](./figures/frame.0100.jpg?raw=true "Stress")|


```
-DBUILD_STRESS=ON 
```

### io 

Reader for the [LAMMPS](https://lammps.sandia.gov/) [file format](https://lammps.sandia.gov/doc/2001/data_format.html)   

```
-DBUILD_IO=ON 
```
## Build instructions

### Dependencies

* [Paraview](https://www.paraview.org/) most be build from source
* [CMake](https://cmake.org/) >= 2.8

### Building the plugins

```
cd filter
mkdir build && cd build
cmake -DParaView_DIR=/calculs/Compile/ParaView-v5.3.0/build \
      -DCMAKE_INSTALL_PREFIX=path/to/install/plugins
make -j
```
### Load the plugins to paraview

```
export PV_PLUGIN_PATH=$PV_PLUGIN_PATH:path/to/install/plugins
```


 

## Publications

1. P. Diehl, M. Bußler, D. Pflüger, S. Frey, T. Ertl, F. Sadlo, and M. A. Schweitzer, Extraction of Fragments and Waves After Impact Damage in Particle-Based Simulations, Springer International Publishing, Cham, 2017, pp. 17–34.
