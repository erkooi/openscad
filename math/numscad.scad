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
// Author: E. Kooistra, May 2022
// Purpose: Support functions from numpy.py in OpenSCAD
// Description:
//   Functions:
//   . num_range(start, step, end)
//   . num_slice(x_arr, start, end)
//   . num_linspace(start, stop, N)
//   . num_sum(x_arr)
//   . num_add_vectors(x_arr, y_arr)
//   . num_add_vector_skalar(x_arr, y)
//   . num_any(b_arr)
//   . num_zero_if_close(x_arr, f_eps=0)
//   . num_isclose(a_arr, b_arr, rtol=1e-05, atol=1e-08)
//   . num_allclose(a_arr, b_arr, rtol=1e-05, atol=1e-08)
//   . num_complex_abs(re_arr, im_arr)
//   . num_complex_angle(re_arr, im_arr, f_eps=0)
//   . num_cos_signal(t_arr, K, phi, dc)
//   . num_sin_signal(t_arr, K, phi, dc)
//
//   Tests:
//   . matrix - matrix multiplication of arbitrary dimensions
//   . matrix - vector multiplication of arbitrary dimensions

//------------------------------------------------------------------------------
// Vector generation
//------------------------------------------------------------------------------

// Generate range vector
function num_range(start, step, end) = [ for (i = [start : step : end]) i ];


// Slice range from vector
function num_slice(x_arr, start, end) = end < start ? [] : [ for (i = [start : end]) x_arr[i] ];


// Make linspace using recursive (_rcs) function and concat(), out of interest
function num_linspace_recursive(start, stop, N) =
    let(step = (stop - start) / N)
    N > 0 ? concat([start], num_linspace(start + step, stop, N-1)) : [];


// Make linspace using list comprehension with for(), default approach
function num_linspace(start, stop, N) =
    let (step = (stop - start) / N)
    num_range(start, step, stop - step);


//------------------------------------------------------------------------------
// Vector operation
// . use list comprehension to determin a vector result per vector element
//   operation
// . use recursion to determine single value result for a vector operation
//------------------------------------------------------------------------------
// Return sum of elements in vector
function num_sum(x_arr) =
    let (N = len(x_arr))
    N == 0 ? 0 :
    N == 1 ? x_arr[0] :
             x_arr[N-1] + num_sum(num_slice(x_arr, 0, N-2));


// Return addition of two vectors
function num_add_vectors(x_arr, y_arr) =
    let (N = len(x_arr),
         n_arr = num_range(0, 1, N-1)
        )
    (len(y_arr) == N)
         ? [ for (n = n_arr) x_arr[n] + y_arr[n] ]
         : false;  // cannot add different length vectors


// Return addition of vector and skalar
function num_add_vector_skalar(x_arr, y) =
    let (N = len(x_arr),
         n_arr = num_range(0, 1, N-1)
        )
    [ for (n = n_arr) x_arr[n] + y ];


// Return boolean vector with not element for each element in vector, similar
// as numpy.num_logical_not()
function num_logical_not(b_arr) =
    let (N = len(b_arr),
         n_arr = num_range(0, 1, N-1)
        )
    [ for (n = n_arr) !b_arr[n] ];


// Return true when any element in vector is true, similar as numpy.any()
function num_any(b_arr) =
    let (N = len(b_arr))
    N == 0 ? false :
    N == 1 ? b_arr[0] :
             b_arr[N-1] || num_any(num_slice(b_arr, 0, N-2));


// Return vector with elements forced to zero when abs(element) < f_eps
function num_zero_if_close(x_arr, f_eps=0) =
    let (N = len(x_arr),
         n_arr = num_range(0, 1, N-1)
        )
    [ for (n = n_arr) (f_eps > 0 && abs(x_arr[n]) < f_eps)
                          ? 0
                          : x_arr[n] ];


// Return boolean vector with element true when two vectors elements are close
// to equal, similar as numpy.isclose(). Uses a_arr as reference for length N
// and compare elements 0:N-1 from b_arr, so len(b_arr) must be >= N.
function num_isclose(a_arr, b_arr, rtol=1e-05, atol=1e-08) =
    let (N = len(a_arr),
         n_arr = num_range(0, 1, N-1),
         // evaluate abs(a - b) <= (atol + rtol * abs(b))
         abs_diff_arr = [ for (n = n_arr) abs(a_arr[n] - b_arr[n]) ],
         tol_diff_arr = [ for (n = n_arr) atol + rtol * abs(a_arr[n]) ]
        )
    // isclose() vector
    [ for (n = n_arr) abs_diff_arr[n] <= tol_diff_arr[n] ];


// Return value true when all elements in two vectors are close to equal,
// similar as numpy.allclose(). Return false if a_arr and b_arr differ in
// length.
function num_allclose(a_arr, b_arr, rtol=1e-05, atol=1e-08) =
    (len(a_arr) == len(b_arr))
        ?   !num_any(num_logical_not(num_isclose(a_arr, b_arr, rtol, atol)))
        :   false;

// Return magnitudes vector of complex real vector and imag vector pair
function num_complex_abs(re_arr, im_arr) =
    let (N = len(re_arr),
         n_arr = num_range(0, 1, N-1)
        )
    [ for (n = n_arr) sqrt(re_arr[n]^2 + im_arr[n]^2) ];


// Return angles vector of complex real vector and imag vector pair
// . If abs(re) < f_eps and abs(im) < f_eps then force angle to 0, like for
//   atan2(0, 0).
function num_complex_angle(re_arr, im_arr, f_eps=0) =
    let (N = len(re_arr),
         n_arr = num_range(0, 1, N-1)
        )
    [ for (n = n_arr) (f_eps > 0 && abs(re_arr[n]) < f_eps && abs(im_arr[n]) < f_eps)
                          ? 0
                          : atan2(im_arr[n], re_arr[n]) ];


//------------------------------------------------------------------------------
// Gonio signals
// . t_arr: time array with N points in [0, 1>
// . K: number of periods in t_arr, for frequency in bin k = K of DFT and rDFT
// . ampl: amplitude
// . phi: angle in degrees
// . dc: dc offset
//------------------------------------------------------------------------------

// Return cos signal using recursive function and concat(), out of interest
function num_cos_signal_rcs(t_arr, K, ampl, phi, dc=0) =
    let (N = len(t_arr),
         t_remainder = num_slice(t_arr, 1, N-1)
        )
    N == 0 ? [] : concat([ampl * cos(360 * t_arr[0] * K + phi) + dc], num_cos_signal_rcs(t_remainder, K, ampl, phi, dc));


// Return cos signal using list comprehension with for(), default approach
function num_cos_signal(t_arr, K, ampl, phi, dc=0) =
    [ for (t = t_arr) ampl * cos(360 * t * K + phi) + dc ];


// Return sin signal using list comprehension with for(), default approach
function num_sin_signal(t_arr, K, ampl, phi, dc=0) =
    [ for (t = t_arr) ampl * sin(360 * t * K + phi) + dc ];


//------------------------------------------------------------------------------
// UNIT TEST
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Test range()
//------------------------------------------------------------------------------
r_incr = num_range(0, 1, 10);
r_decr = num_range(0, -3, -10);
r_empty = num_range(0, -3, 10);
echo();
echo(">>> Test range()");
echo(r_incr = r_incr);
echo(r_decr = r_decr);
echo(r_empty = r_empty);

assert (r_incr == [0, 1, 2, 3, 4, 5, 6, 7 , 8, 9, 10], "Wrong r_incr.");
assert (r_decr == [0, -3, -6, -9], "Wrong r_decr");
assert (r_empty == [], "Wrong r_empty != []");


//------------------------------------------------------------------------------
// Test slice()
//------------------------------------------------------------------------------
N_slice = 10;  // fixed 10 for assert tests
i_arr = num_range(0, 1, N_slice-1);
i0 = num_slice(i_arr, 0, 0);
i1 = num_slice(i_arr, 1, 1);
i_empty = num_slice(i_arr, 0, -1);
i_slice = num_slice(i_arr, N_slice-2, N_slice-1);
i_slice0 = i_slice[0];   // sliced result starts at index 0

echo();
echo(">>> Test slice()");
echo(i_arr = i_arr);
echo(i0 = i0);  // = i_arr[0]
echo(i1 = i1);  // = i_arr[1]
echo(i_empty = i_empty);  // empty list for end < start
echo(i_slice = i_slice);  // = i_arr[slice]
echo(i_slice0 = i_slice0);

assert (i_empty == [], "Wrong i_empty");
assert (i_slice == [8, 9], "Wrong i_slice");
assert (i_slice0 == 8, "Wrong i_slice0");


//------------------------------------------------------------------------------
// Test linspace()
//------------------------------------------------------------------------------
// . Use N = 8 is power of 2 to avoid float quantization differences in == test.
ls_rcs = num_linspace_recursive(0, 1, 8);
ls_arr = num_linspace(0, 1, 8);
echo();
echo(">>> Test linspace()");
echo(ls_rcs = ls_rcs);
echo(ls_arr = ls_arr);

assert(ls_rcs == ls_arr, "ls_rcs != ls_arr");


//------------------------------------------------------------------------------
// Test sum()
//------------------------------------------------------------------------------
empty_sum = num_sum([]);

n_arr = num_range(0, 1, 10);

n_sum = num_sum(n_arr);

echo();
echo(">>> Test sum");
echo(empty_sum = empty_sum);
echo(n_arr = n_arr);
echo(n_sum = n_sum);

assert (empty_sum == 0, "Wrong empty_sum");
assert (n_sum == 55, "Wrong n_sum");

//------------------------------------------------------------------------------
// Test add_vectors(), add_vector_skalar()
//------------------------------------------------------------------------------
x_arr = [0, 1, 2, 3];
y_arr = [4, 5, 6, 7];
e_arr = [3, 2];

echo();
echo(">>> Test add_vectors(), add_vector_skalar()");
assert (num_add_vector_skalar(x_arr, 1) == [1, 2, 3, 4], "Wrong add vector + skalar");
assert (num_add_vectors(x_arr, y_arr) == [4, 6, 8, 10], "Wrong add two vectors");
assert (num_add_vectors(x_arr, e_arr) == false, "Wrong add two vectors");


//------------------------------------------------------------------------------
// Test logical_not()
// . !arr yields a boolean skalar, and ![false, false, false] yields false,
//   therefore use ! for skalar and define num_logical_not() for vector.
//------------------------------------------------------------------------------
all_true_arr = [true, true, true];
all_false_arr = [false, false, false];
some_true_arr = [false, true, false];
some_false_arr = [true, false, true];  // = not some_true_arr
not_some_true_arr = num_logical_not(some_true_arr);
not_all_false_arr = num_logical_not(all_false_arr);

echo();
echo(">>> Test logical_not()");
echo(some_true_arr = some_true_arr)
echo(not_some_true_arr = not_some_true_arr)
echo(all_false_arr = all_false_arr)
echo(not_all_false_arr = not_all_false_arr)

assert (not_some_true_arr == some_false_arr, "Wrong num_logical_not");
assert (not_all_false_arr == all_true_arr, "Wrong num_logical_not");


//------------------------------------------------------------------------------
// Test any()
//------------------------------------------------------------------------------
echo();
echo(">>> Test any()");
assert (num_any(some_true_arr) == true, "Wrong num_any for one or more true array");
assert (num_any(all_false_arr) == false, "Wrong num_any for all false array");


//------------------------------------------------------------------------------
// Test zero_if_close()
//------------------------------------------------------------------------------

almost_zero_arr = [0.1, 1e-5, -1e-3, 0];

echo();
echo(">>> Test zero_if_close()");
assert (num_zero_if_close(almost_zero_arr) == almost_zero_arr, "Wrong num_zero_if_close");
assert (num_zero_if_close(almost_zero_arr, f_eps=1e-1) == [0.1, 0, 0, 0], "Wrong num_zero_if_close");
assert (num_zero_if_close(almost_zero_arr, f_eps=1e-4) == [0.1, 0, -1e-3, 0], "Wrong num_zero_if_close");


//------------------------------------------------------------------------------
// Test zero_if_close, isclose(), allclose()
//------------------------------------------------------------------------------
a_arr = [1.002, -3.001, 5.02];
b_arr = [1.0, -3.0, 5.0];
c_arr = [1.0, -3.0, 5.0, 6.1];   // different lenght than b_arr

echo();
echo(">>> Test isclose(), allclose()");
echo(a_arr = a_arr);
echo(b_arr = b_arr);
echo(c_arr = c_arr);
echo("num_isclose rtol=0.01:", num_isclose(a_arr, b_arr, rtol=0.01));
echo("num_isclose rtol=0.001:", num_isclose(a_arr, b_arr, rtol=0.001));
echo("num_isclose atol=0.1:", num_isclose(a_arr, b_arr, atol=0.1));
echo("num_isclose atol=0.01:", num_isclose(a_arr, b_arr, atol=0.01));
echo("num_isclose:", num_isclose(b_arr, c_arr));

// . tolerances
assert (num_allclose(a_arr, b_arr, rtol=0.01) == true, "Wrong num_allclose for rtol = 0.01");
assert (num_allclose(a_arr, b_arr, rtol=0.001) == false, "Wrong num_allclose for rtol = 0.001");
assert (num_allclose(a_arr, b_arr, atol=0.1) == true, "Wrong num_allclose for atol = 0.1");
assert (num_allclose(a_arr, b_arr, atol=0.01) == false, "Wrong num_allclose for atol = 0.01");

// . different lengths
assert (num_isclose(a_arr, c_arr, atol=0.01) == [true, true, false], "Wrong num_isclose, for different length");
assert (num_isclose(b_arr, c_arr) == [true, true, true], "Wrong num_isclose, for different length");
assert (num_allclose(b_arr, c_arr) == false, "Wrong num_allclose for different length");


//------------------------------------------------------------------------------
// Test cos_signal()
//------------------------------------------------------------------------------
N_time = 32;
t_arr = num_linspace(0, 1, N_time);
K = 1;  // one period in t_arr
ampl = 2;
phi = 30;
dc = 0.1;

cos_rcs = num_cos_signal_rcs(t_arr, K, ampl, phi, dc);
cos_arr = num_cos_signal(t_arr, K, ampl, phi, dc);
sin90_arr = num_sin_signal(t_arr, K, ampl, phi + 90, dc);  // sin(phi + 90) = cos(phi)
cos_empty_rcs = num_cos_signal_rcs([], K, ampl, phi, dc);
cos_empty_arr = num_cos_signal([], K, ampl, phi, dc);

echo();
echo(">>> Test cos_signal()");
echo(cos_rcs = cos_rcs);
echo(cos_arr = cos_arr);
echo(sin90_arr = sin90_arr);
echo(cos_empty_rcs = cos_empty_rcs);
echo(cos_empty_arr = cos_empty_arr);

assert (cos_rcs == cos_arr, "Wrong cos_rcs != cos_arr");
assert (cos_arr == sin90_arr, "Wrong cos_arr != sin90_arr");
assert (cos_empty_rcs == [], "Wrong cos_empty_rcs");
assert (cos_empty_arr == [], "Wrong cos_empty_arr");


//------------------------------------------------------------------------------
// Test num_complex_*()
//------------------------------------------------------------------------------
test_re_arr = [0, 1, 1];
test_im_arr = [1, 1, 0];

test_abs = num_complex_abs(test_re_arr, test_im_arr);
test_angle = num_complex_angle(test_re_arr, test_im_arr);

echo();
echo(">>> Test complex abs(), angle()");
echo(test_abs = test_abs);
echo(test_angle = test_angle);

assert (test_abs == [1, sqrt(2), 1], "Wrong test_abs");
assert (test_angle == [90, 45, 0], "Wrong test_angle");


//------------------------------------------------------------------------------
// Try:
// . Matrix - matrix multiplication
// . Matrix - vector multiplication
// . Vector - vector multiplication is not possible with standard OpenSCAD *
//            operator
// . Matrix - skalar multiplication
// . Vector - skalar multiplication
// . Vector - vector, skalar addition is not possible with standard OpenSCAD +
//            operator, therefore use num_add_vectors()
//------------------------------------------------------------------------------

M2 = [ [1, 2],
       [3, 4] ];
I2 = [ [1, 0],
       [0, 1] ];
M23 = [ [1, 2, 3],
        [4, 5, 6] ];
M3 = [ [1, 2, 3],
       [4, 5, 6],
       [7, 8, 9] ];
I3 = [ [1, 0, 0],
       [0, 1, 0],
       [0, 0, 1] ];
x2 = [1, 0];
x3 = [0, 1, 0];
echo();
echo(">>> Try Matrix * Matrix");
echo("M2 * I2", M2 * I2);
echo();
echo(">>> Try Matrix * Vector");
echo("M2 * x2", M2 * x2);
echo("M23 * x3", M23 * x3);
echo("M3 * x3", M3 * x3);
echo(">>> Try skalar *");
echo("M3 * 3", M3 * 3);
echo("x3 * 3", x3 * 3);

M = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];
echo();
echo(">>> Try Matrix row column order: M[row][col] as with python");
echo(M = M);
echo("M[0][1] = ", M[0][1]);
echo("M[1][0] = ", M[1][0]);


//------------------------------------------------------------------------------
// All asserts when OK when program gets to here
//------------------------------------------------------------------------------
echo("PASSED");