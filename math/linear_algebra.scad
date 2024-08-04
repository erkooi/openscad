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
// Purpose: Library of Linear Algebra functions
// Description:
//   Contains defintions and functions for:
//    . coordinates system x, y, z
//    . rotation matrices that are equivalent to OpenSCAD rotate(),
//    . translation matrices that are equivalent to OpenSCAD translate(),
//    . angles using atan2, but with atan2(0, 0) = 0 returned as undef to show
//      that the angle is undefined (can be any value), because the point lies
//      on the rotation axis (or within some almost zero distance, due to
//      floating point accuracy).
// Remark:
// . uses prefix 'm' for OpenSCAD module
// . uses prefix 'f' for OpenSCAD function

//------------------------------------------------------------------------------
// Coordinates
//------------------------------------------------------------------------------
//
// Right-handed coordinates X, Y, Z:
// . right hand fingers point from X to Y points tumb to Z
// . screw from X to Y goes to Z
//
//     z
//     |          . angleYZ
//     |--- y     . angleZX
//    /           . angleXY
//   x
//
// Euler angles
// - https://en.wikipedia.org/wiki/Euler_angles
// - positive rotation about:
//   . X-axis from Y to Z = phi, is roll, aileron, banking
//   . Y-axis from Z to X = theta, is pitch, elevator, elevation
//   . Z-axis from X to Y = psi, is yaw, rudder, heading
//
// Spherical coordinates P(r, phi, theta)
// - https://en.wikipedia.org/wiki/Spherical_coordinate_system
//   . r is radial distance from origin to P
//   . phi is angle from X to projection of P in XY plane
//   . theta is angle from Z to radial line through P

//------------------------------------------------------------------------------
// Matrices
//------------------------------------------------------------------------------

// Translate
module mTrans(x, y, z) {
   multmatrix(m = fTrans(x, y, z)) children();
}

// Rotate
module mRot_yz(phi) {
    multmatrix(m = fRot_yz(phi)) children();  // positive angle from +y to +z around x-axis
}

module mRot_zx(theta) {
    multmatrix(m = fRot_zx(theta)) children();  // positive angle from +z to +x around y-axis
}

module mRot_xy(psi) {
    multmatrix(m = fRot_xy(psi)) children();  // positive angle from +x to +y around z-axis
}

// Translation matrices
function fTrans(x, y, z) = [ [ 1, 0, 0, x],
                             [ 0, 1, 0, y],
                             [ 0, 0, 1, z],
                             [ 0, 0, 0, 1] ];

// Rotation axis matrices
function fRot_yz(phi) = [ [ 1,        0,         0, 0],
                          [ 0, cos(phi), -sin(phi), 0],
                          [ 0, sin(phi),  cos(phi), 0],
                          [ 0,        0,         0, 1] ];  // positive angle from +y to +z around x-axis

function fRot_zx(theta) = [ [ cos(theta), 0,  sin(theta), 0],
                            [          0, 1,           0, 0],
                            [-sin(theta), 0,  cos(theta), 0],
                            [          0, 0,           0, 1] ];  // positive angle from +z to +x around y-axis

function fRot_xy(psi) = [ [ cos(psi), -sin(psi), 0, 0],
                          [ sin(psi),  cos(psi), 0, 0],
                          [        0,         0, 1, 0],
                          [        0,         0, 0, 1] ];  // positive angle from +x to +y around z-axis

//------------------------------------------------------------------------------
// Angles for vector p = [x, y, z, w]
//------------------------------------------------------------------------------

// Rotation axis angles
// . atan2(b, a) = arctan(b / a) in range [-180, 180] degrees
// . atan2() can become -180 or +180 when perpendicular side is -0 or +0
// . atan2(0, 0) = 0, so undefined angle forced to 0, like arctan2() in numpy
// . xyz coordinates, positive angles from +x to +y to +z to +x:
//   . Ryz: positive angle from +y to +z around x-axis, angle = atan2(z, y)
//   . Rzx: positive angle from +z to +x around y-axis, angle = atan2(x, z)
//   . Rxy: positive angle from +x to +y around z-axis, angle = atan2(y, x)
//   . Rxr: positive angle is from +x to r = sqrt(y^2 + z^2), so off axis angle = atan2(r, x)
//   . Ryr: positive angle is from +y to r = sqrt(z^2 + x^2), so off axis angle = atan2(r, y)
//   . Rzr: positive angle is from +z to r = sqrt(x^2 + y^2), so off axis angle = atan2(r, z)

// Return atan2(0, 0) = 0
function fAngleYZ(p) = atan2(p.z, p.y);  // positive angle from +y to +z around x-axis, roll
function fAngleZX(p) = atan2(p.x, p.z);  // positive angle from +z to +x around y-axis, pitch
function fAngleXY(p) = atan2(p.y, p.x);  // positive angle from +x to +y around z-axis, jaw
function fAngleXR(p) = let (r = sqrt(p.y^2 + p.z^2)) atan2(r, p.x);  // positive angle from +x to r off axis
function fAngleYR(p) = let (r = sqrt(p.z^2 + p.x^2)) atan2(r, p.y);  // positive angle from +y to r off axis
function fAngleZR(p) = let (r = sqrt(p.x^2 + p.y^2)) atan2(r, p.z);  // positive angle from +z to r off axis

function toAngle360(angle) = let (a = angle % 360) a < 0 ? a + 360 : a;  // range [0,360>
function toAngle180(angle) = let (a = toAngle360(angle)) a > 180 ? a - 360 : a;  // range <-180, 180]


//------------------------------------------------------------------------------
// UNIT TEST
//------------------------------------------------------------------------------

p = [2, 1, 1, 1];
echo();
echo("Angles:");
echo("fAngleYZ(p) = ", fAngleYZ(p));
echo("fAngleZX(p) = ", fAngleZX(p));
echo("fAngleXY(p) = ", fAngleXY(p));
echo("fAngleXR(p) = ", fAngleXR(p));
echo("fAngleYR(p) = ", fAngleYR(p));
echo("fAngleZR(p) = ", fAngleZR(p));

// Openscad uses doubles but echo and str round to 5 digits
echo();
p1 = [0, -1,  1e-5, 0]; echo(str("fAngleYZ(", p1, ") = ", fAngleYZ(p1)));
p2 = [0, -1, -1e-5, 0]; echo(str("fAngleYZ(", p2, ") = ", fAngleYZ(p2)));
p3 = [0, -1,  1e-6, 0]; echo(str("fAngleYZ(", p3, ") = ", fAngleYZ(p3)));
p4 = [0, -1, -1e-6, 0]; echo(str("fAngleYZ(", p4, ") = ", fAngleYZ(p4)));

// x > 0, x % 360 yields 0 <= value < 360
// x < 0, x % 360 yields -360 < value <= 0
echo();
echo("Modulo:");
echo(" 270 % 360 ",  270 % 360);
echo(" 360 % 360 ",  360 % 360);
echo(" 450 % 360 ",  450 % 360);
echo("-270 % 360 ", -270 % 360);
echo("-360 % 360 ", -360 % 360);
echo("-450 % 360 ", -450 % 360);

// atan2() yields -180 <= angle <= 180
echo();
echo("atan2():");
echo("atan2( 0,  0) = ", atan2( 0,  0));
echo("atan2( 0,  1) = ", atan2( 0,  1));
echo("atan2( 1,  0) = ", atan2( 1,  0));
echo("atan2( 0, -1) = ", atan2( 0, -1));
echo("atan2(-1,  0) = ", atan2(-1,  0));

// fAngleYZ() yields -180 < angle <= 180 and undef when Y, Z = 0, 0
//                 Y  Z
echo();
echo("fAngleYZ([0, 0, 0, 0]) = ", fAngleYZ([0, 0, 0, 0]));
echo("fAngleYZ([0, 1, 0, 0]) = ", fAngleYZ([0, 1, 0, 0]));
echo("fAngleYZ([0, 0, 1, 0]) = ", fAngleYZ([0, 0, 1, 0]));
echo("fAngleYZ([0,-1, 0, 0]) = ", fAngleYZ([0,-1, 0, 0]));
echo("fAngleYZ([0, 0,-1, 0]) = ", fAngleYZ([0, 0,-1, 0]));

// toAngle360(angle) yields 0 <= value < 360, also for angle < 0
echo();
echo("toAngle360():");
echo("   0 = ", toAngle360(   0));
echo(" -90 = ", toAngle360( -90));
echo("-180 = ", toAngle360(-180));
echo("-270 = ", toAngle360(-270));
echo("-360 = ", toAngle360(-360));
echo("-450 = ", toAngle360(-450));
echo("   0 = ", toAngle360(   0));
echo("  90 = ", toAngle360(  90));
echo(" 180 = ", toAngle360( 180));
echo(" 270 = ", toAngle360( 270));
echo(" 360 = ", toAngle360( 360));
echo(" 450 = ", toAngle360( 450));

// toAngle180(angle) yields -180 < value <= 180, like atan2()
echo();
echo("toAngle180():");
echo("   0 = ", toAngle180(   0));
echo(" -90 = ", toAngle180( -90));
echo("-180 = ", toAngle180(-180));
echo("-270 = ", toAngle180(-270));
echo("-360 = ", toAngle180(-360));
echo("-450 = ", toAngle180(-450));
echo("   0 = ", toAngle180(   0));
echo("  90 = ", toAngle180(  90));
echo(" 180 = ", toAngle180( 180));
echo(" 270 = ", toAngle180( 270));
echo(" 360 = ", toAngle180( 360));
echo(" 450 = ", toAngle180( 450));
