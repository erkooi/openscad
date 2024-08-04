# Libraries for OpenSCAD

## 1. Mathematical operations
openscad/math/
* common.scad : Common constants, functions
* linear_algebra.scad : Linear Algebra functions for matrices, angles
* numscad.scad : Support functions from Pyhton numpy in OpenSCAD
* dft.scad : Calculate DFT for real input using matrix multiplication
* triangles.scad : Create a triangle shape marker module object

openscad/math/python/
* libraries/linear_algebra.py : Linear Algebra functions in Python for comparison with linear_algebra.scad
* try_dft.ipynb, try_dft.html : Calculate DFT for a real input using operations that can also be applied in OpenSCAD dft.scad

## 2. Shapes
openscad/shapes
* triangles.scad : Create a triangle shape marker