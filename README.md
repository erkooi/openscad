# Libraries for OpenSCAD

## 1. Mathematical operations
openscad/math/
* math_constants.scad : Math constants to include
* linear_algebra.scad : Linear Algebra functions for matrices, angles
* numscad.scad : Support functions from Python numpy in OpenSCAD
* dft.scad : Calculate DFT for real input using matrix multiplication
* triangles.scad : Create a triangle shape marker module object

openscad/math/python/
* libraries/linear_algebra.py : Linear Algebra functions in Python for comparison with linear_algebra.scad
* try_dft.ipynb, try_dft.html : Calculate DFT for a real input using operations that can also be applied in OpenSCAD dft.scad

## 2. Shapes
openscad/shapes
* shape_constants.scad : Shape constants to include
* triangles.scad : Create a triangle shape marker