################################################################################
# Copyright 2022 E. Kooistra
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
################################################################################
#
# Author: E. Kooistra, March 2022
# Purpose: Library of Linear Algebra functions
# Description:
# . Corresponds to linear_algebra.scad
# Remark:
# . np.arctan2(b, a) = arctan(b / a) in range [-pi, pi]
# . np.arctan2() can become -pi or +pi when perpendicular side is -0 or +0

import numpy as np

################################################################################
# Constants
################################################################################

# from common.scad
f_eps = 1e-10  # Floating point value considered equal to 0


################################################################################
# Matrices
################################################################################

# Translation matrices
def Trans(x, y, z):
    """Translate by [x, y, z]."""
    return np.array([[1, 0, 0, x],
                     [0, 1, 0, y],
                     [0, 0, 1, z],
                     [0, 0, 0, 1]])


# Rotation axis matrices
def Rot_yz(angle):
    """Rotate around x-axis from y to z."""
    a = np.radians(angle)
    co, si = np.cos(a), np.sin(a)
    return np.array([[1, 0, 0, 0],
                     [0, co, -si, 0],
                     [0, si, co, 0],
                     [0, 0, 0, 1]])


def Rot_zx(angle):
    """Rotate around y-axis from z to x."""
    a = np.radians(angle)
    co, si = np.cos(a), np.sin(a)
    return np.array([[co, 0, si, 0],
                     [0, 1, 0, 0],
                     [-si, 0, co, 0],
                     [0, 0, 0, 1]])


def Rot_xy(angle):
    """Rotate around z-axis from x to y."""
    a = np.radians(angle)
    co, si = np.cos(a), np.sin(a)
    return np.array([[co, -si, 0, 0],
                     [si, co, 0, 0],
                     [0, 0, 1, 0],
                     [0, 0, 0, 1]])


################################################################################
# Angles for vector p = [x, y, z, w]
################################################################################

def fAngleYZ(p):
    """Roll angle from y to z around x-axis."""
    y = p[1]
    z = p[2]
    if np.abs(y) < f_eps and np.abs(z) < f_eps:
        return np.nan
    else:
        return np.degrees(np.arctan2(z, y))


def fAngleZX(p):
    """Pitch angle from z to x around y-axis."""
    x = p[0]
    z = p[2]
    if np.abs(x) < f_eps and np.abs(z) < f_eps:
        return np.nan
    else:
        return np.degrees(np.arctan2(x, z))


def fAngleXY(p):
    """Jaw angle from x to y around z-axis."""
    x = p[0]
    y = p[1]
    if np.abs(x) < f_eps and np.abs(y) < f_eps:
        return np.nan
    else:
        return np.degrees(np.arctan2(y, x))


def fAngleArrYZ(angleArr):
    """Roll angles from y to z around x-axis."""
    return np.array([fAngleYZ(x) for x in angleArr])


def fAngleArrZX(angleArr):
    """Pitch angles from z to x around y-axis."""
    return np.array([fAngleZX(x) for x in angleArr])


def fAngleArrXY(angleArr):
    """Jaw angles from x to y around z-axis."""
    return np.array([fAngleXY(x) for x in angleArr])


def fAngleXR(p):
    """Angle from x to r around origin."""
    x = p[0]
    y = p[1]
    z = p[2]
    r = np.sqrt(np.power(z, 2) + np.power(y, 2))
    return np.degrees(np.arctan2(r, x))


def fAngleYR(p):
    """Angle from y to r around origin."""
    x = p[0]
    y = p[1]
    z = p[2]
    r = np.sqrt(np.power(z, 2) + np.power(x, 2))
    return np.degrees(np.arctan2(r, y))


def fAngleZR(p):
    """Angle from z to r around origin."""
    x = p[0]
    y = p[1]
    z = p[2]
    r = np.sqrt(np.power(x, 2) + np.power(y, 2))
    return np.degrees(np.arctan2(r, z))


def fAngleArrXR(angleArr):
    """Angles from x to r around x-axis."""
    return np.array([fAngleXR(x) for x in angleArr])


def fAngleArrYR(angleArr):
    """Angles from y to r around y-axis."""
    return np.array([fAngleYR(x) for x in angleArr])


def fAngleArrZR(angleArr):
    """Angles from z to r around z-axis."""
    return np.array([fAngleZR(x) for x in angleArr])


def toAngle360(angle):
    """Map angle to range [0, 360> degrees."""
    # Check input nan, to avoid RuntimeWarning: invalid value encountered in
    # double_scalars with nan % 360, that can occur with toAngleArr360() in
    # jupyter notebook.
    if np.isnan(angle):
        return np.nan
    a = angle % 360
    if a < 0:
        return a + 360
    else:
        return a


def toAngle180(angle):
    """Map angle to range <-180, 180] degrees."""
    a = toAngle360(angle)
    if a > 180:
        return a - 360
    else:
        return a


def toAngleArr360(angleArr):
    """Map array of angles to range [0, 360> degrees."""
    return np.array([toAngle360(x) for x in angleArr])


def toAngleArr180(angleArr):
    """Map array of angles to range <-180, 180] degrees."""
    return np.array([toAngle180(x) for x in angleArr])


if __name__ == '__main__':

    # Use @ to apply sequency of matrices
    a = np.array([0, 0, 0, 1])
    b = Rot_zx(90) @ Rot_yz(90) @ Rot_xy(90) @ Trans(15, 0, 0) @ a
    print('Rot_zx(90) @ Rot_yz(90) @ Rot_xy(90) @ Trans(15, 0, 0) @', a, '=', b)
    print()

    # Use same values as in linear_algebra.scad to check that they yield the same result
    p = [2, 1, 1, 1]
    print('p = ', p)
    print('fAngleYZ(p) = ', fAngleYZ(p))
    print('fAngleZX(p) = ', fAngleZX(p))
    print('fAngleXY(p) = ', fAngleXY(p))
    print('fAngleXR(p) = ', fAngleXR(p))
    print('fAngleYR(p) = ', fAngleYR(p))
    print('fAngleZR(p) = ', fAngleZR(p))

    # arctan2() for YZ can become -180 or 180 when Z is -0 or +0
    print()
    print('fAngleYZ([0, -1,  1e-15, 0]) = ', fAngleYZ([0, -1,  1e-15, 0]))
    print('fAngleYZ([0, -1, -1e-15, 0]) = ', fAngleYZ([0, -1, -1e-15, 0]))
    print('fAngleYZ([0, -1,  1e-16, 0]) = ', fAngleYZ([0, -1,  1e-16, 0]))
    print('fAngleYZ([0, -1, -1e-16, 0]) = ', fAngleYZ([0, -1, -1e-16, 0]))

    print()
    print('toAngle360( 360) = ', toAngle360(360))
    print('toAngle360(-1)   = ', toAngle360(-1))
    print('toAngle360(-181) = ', toAngle360(-181))
    print('toAngle360(-180) = ', toAngle360(-180))
    print('toAngle360( 180) = ', toAngle360(180))
