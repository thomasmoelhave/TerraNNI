TerraNNI - Natural Neighbor Interpolation on 2D and 3D point clouds
====================================================================

TerraNNI interpolates large 2D and 3D pointclouds using a variant of
natural neighbor interpolation. The original natural neighbor
interpolation scheme[1] is approximated by discretizing the Voronoi
diagram and by limiting the region of influence of the sites. More
details in [2] for the 2D version and [3] for the 3D version.

Authors
=======

TerraNNI is developed by:

Pankaj K. Agarwal, Duke University - http://www.cs.duke.edu/~pankaj

Alex Beutel, Carnegie Mellon University - http://alexbeutel.com/

Thomas Mølhave, Duke University - http://www.cs.duke.edu/~thomasm/


Dependencies
============

TerraNNI uses the following libraries

  * libLAS 1.6, available at http://liblas.org/ 
  * TPIE, available at https://github.com/thomasmoelhave/tpie (tpie is
    included through a git submodule, so there is no need to install
    it independently).
  * Boost, available at http://boost.org
  * CUDA 5, available at http://nvidia.com


Compilation
===============

TerraNNI uses CMake for its build system.

```bash
$ git clone git://github.com/thomasmoelhave/TerraNNI.git
$ git submodule init
$ git submodule update
$ mkdir build
$ cd build
$ ccmake ../TerraNNI #( alternatively, use cmake-gui ../TerraNNI )
$ make
```


Caveats
-------

A couple of the CMake variables are vital for the compilation to
succeede.  **CUDA_SDK_ROOT_DIR** must be set to the *C* directory of the
NVIDIA GPU Computing SDK.

If you are on a system running GCC 4.5 or newer, you may have to force
the cuda compiler (nvcc) to use an older version of gcc. This can be
done by setting the CUDA_HOST_COMPILER cmake variable to g++-4.4.



Running
=======

TerraNNI needs a directory for storing temporary files, this is set
with the --tmp-path parameter.  The following command works on the
point cloud build from the las files available at --tilepath. It
compites a 2000 by 2000 grid with a cell size of 1. It starts at year
2002 and stops at year 2004, outputting a grid for every half year in
the mean time. The spatial origin of the output grids is (10,10).

**Note**: TerraNNI interprets the "intensity" value for LAS points as
  the time stamp for that point. Thus, in the example below, the
  intensity values are presumably somewhere in the [2000,2006] range.

```bash
$ ./terranni --tilepath=<path to las files> --output=out_directory --cell-size=1.0 --site-radius=5.0 --wp-path=/var/tmp/wp --time-start=2002 --time-length=3 --time-radius=2 --tmp-path=/var/tmp  --time-step=0.5 --origin-x=10 --origin-y=10 --grid-cols=2000 --grid-rows=200
```

The --site-radius and --time-radious parameters set the *region of
influence* in space and time respectively. When interpolation at a
location, points outside its region of influence is not considered,
consult[2,3] for more details.

Creating traditional 2D grids
------------------------------

The program terranni-2d is produced during compilation. It is equivalent to the original terranni program except that all options concerning time has been disabled, it is therefore somewhat simpler to use for creating standard DEMs.

References
===========

[1] Sibson, R. (1981). "A brief description of natural neighbor interpolation (Chapter 2)". In V. Barnett. Interpreting Multivariate Data. Chichester: John Wiley. pp. 21–36   
[2] Alex Beutel, Thomas Mølhave, Pankaj K. Agarwal (2010) Natural neighbor interpolation based grid DEM construction using a GPU In GIS '10: Proceedings of the 18th ACM SIGSPATIAL International Symposium on Advances in Geographic Information Systems.   
[3] Alex Beutel, Thomas Mølhave, Pankaj K. Agarwal, Arnold P. Boedihardjo, James A. Shine (2011). "TerraNNI: Natural Neighbor Interpolation on a 3D Grid Using a GPU". In GIS '11 Proceedings of the 19th ACM SIGSPATIAL International Symposium on Advances in Geographic Information Systems, 2011.
