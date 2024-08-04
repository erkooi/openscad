//------------------------------------------------------------------------------
// Copyright 2022 E. Kooistra
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//------------------------------------------------------------------------------
//
// Author: E. Kooistra, March 2022
// Purpose: Library to create a triangle shape marker.
// Description:
//  . The marker can be place on an element to provide an orientation marker.

//------------------------------------------------------------------------------
// TriangleCylinder
//------------------------------------------------------------------------------

module TriangleCylinder(R_bottom, R_top, h) {
    // Create triangle volume with 3 equal edges using a cylinder with 3 facets:
    //
    //     cylinder($fn = 3, h=h, r1=R_bottom, r2=R_top, center=false);
    //
    // This is equivalent to polyhedron() of points and faces. Centred with
    // radius R and one corner on X-axis at [R_bottom, 0, 0].

    TrianglePoints = [
        [   R_bottom,                   0, 0],  //0
        [-R_bottom/2,  R_bottom*sqrt(3)/2, 0],  //1
        [-R_bottom/2, -R_bottom*sqrt(3)/2, 0],  //2
        [      R_top,                   0, h],  //3
        [   -R_top/2,     R_top*sqrt(3)/2, h],  //4
        [   -R_top/2,    -R_top*sqrt(3)/2, h]]; //5

    TriangleFaces = [
        [0,1,2],   // bottom
        [0,2,5,3], // side
        [2,1,4,5], // side
        [0,3,4,1], // side
        [3,5,4]];  // top

    polyhedron(TrianglePoints, TriangleFaces);
}

//------------------------------------------------------------------------------
// TriangleMarker
//------------------------------------------------------------------------------

module TriangleMarker(Lx_bottom, Wy_bottom, Wy_top, Hz) {
    // Create triangle with base edge of width Wy_bottom at the Y-axis and two equal
    // length edges that meet in a corner at [Lx_bottom, 0, 0]
    //
    // Translate TriangleCylinder to have back edge on Y-axis and then
    // resize to have X corner position to [Lx_bottom, 0, 0]

    resize([Lx_bottom, Wy_bottom, Hz]) translate([Wy_bottom/2, 0, 0]) TriangleCylinder(Wy_bottom, Wy_top, Hz);
}


//------------------------------------------------------------------------------
// UNIT TEST
//------------------------------------------------------------------------------

TriangleMarker(6, 3, 1, 1);
