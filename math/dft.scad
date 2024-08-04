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
// Purpose: Calculate Discrete Fourier Transform (DFT) for real input, to be
//   able to determine the frequency components of a signal (or curve) directly
//   in OpenSCAD.
// Description:
// - See try_dft.ipynb
// - DFT : for complex input
//   . matrix[row][col], so for k (frequency) in rows, for n (time) in columns
//
//                     N-1       -j 2pi/N kn
//   . X(k) = DFT(x) = sum x(n) e
//                     n=0
//
//                     N-1
//                   = sum x(n) [ cos(2pi/N kn) - j sin(2pi/N kn) ]
//                     n=0
//
//   . X =  W x
//
//               -j2pi/N
//          w = e
//
//          wn = w**n
//
//              | 1  1      1       1       ... 1           |
//              | 1  w      w2      w3      ... w(N-1)      |
//              | 1  w2     w4      w6      ... w2(N-1)     |
//              | 1  w3     w6      w9      ... w3(N-1)     |
//          W = | ...                                       |
//              | 1  w(N-1) w2(N-1) w3(N-1) ... w(N-1)(N-1) |
//
//          W = W_re + j * W_im
//
//   . Note: w**(nk) = w**(nk % N), because w**N = exp(-j2pi) = 1
//
// - rDFT : DFT for real input
//
// - Matrix format [complex][K][N], where
//   . complex is [re, im]
//   . K is frequency, k = 0 : K-1, with K = floor(N / 2) + 1
//   . N is time sample, n = 0 : N - 1
//   so:
//       [0] = real matrix[row][col]
//       [1] = imag matrix[row][col]
//       for k (frequency) in rows, for n (time) in columns
//
// - Note:
//   Implements the rDFT matrix multiplication, so the number of operations
//   increases with N^2, and not with N log(N) as with the FFT (Fast Fourier
//   Transform). This can make simulations and animations slow.

use <numscad.scad>;

// rDFT matrix W, for real input
function rdft_W(N) = let(u = 360 / N,
                         K = floor(N / 2) + 1,
                         re = [ for (k = [0 : 1 : K-1]) [ for (n = [0 : 1 : N-1]) cos(-u * ((k * n) % N)) ] ],
                         im = [ for (k = [0 : 1 : K-1]) [ for (n = [0 : 1 : N-1]) sin(-u * ((k * n) % N)) ] ])
                     [re, im];

// rDFT(x) = X = rdft_W * x, for real input
function rdft(x_arr) =
    let (N = len(x_arr),
         W = rdft_W(N),
         re = W[0] * x_arr,
         im = W[1] * x_arr)
    [re, im];


// DFT scale
function dft_scale_dc_harmonics(X, N) =
    // . scale k = 0 DC by 1 / N
    // . scale k > 0 harmonics by 2 / N to get their amplitude
    let(K = len(X))
    concat(X[0], num_slice(X, 1, K - 1) * 2) / N;


// DFT shift, only used with complex input DFT
function dft_fftshift(X) =
    let (N = len(X), K = N % 2 == 0 ? N / 2 : (N-1) / 2)
    N == 0 ? [] :
    N == 1 ? X :
             concat([ for (i = [K : N-1]) X[i] ], [ for (i = [0 : K-1]) X[i] ]);


//------------------------------------------------------------------------------
// UNIT TEST
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Test rDFT matrix W
//------------------------------------------------------------------------------
N_test = 4;
Wtest = rdft_W(N_test);

echo();
echo(">>> Test rDFT matrix W");
echo(N_test = N_test);
echo("Wtest real =", Wtest[0]);
echo("Wtest imag =", Wtest[1]);

if (N_test == 4) {
    assert (Wtest[0] == [[1,  1,  1,  1],
                         [1,  0, -1,  0],
                         [1, -1,  1, -1]], "Wrong Wtest real");

    assert (Wtest[1] == [[0,  0,  0,  0],
                         [0, -1,  0,  1],
                         [0,  0,  0,  0]], "Wrong Wtest imag");
}


//------------------------------------------------------------------------------
// Test dft_scale_dc_harmonics
//------------------------------------------------------------------------------
echo();
echo(">>> Test dft_scale_dc_harmonics");
assert ([1, 2, 2, 2] == dft_scale_dc_harmonics([4, 4, 4, 4], 4), "Wrong dft_scale_dc_harmonics");


//------------------------------------------------------------------------------
// Test dft_fftshift
//------------------------------------------------------------------------------
echo();
echo(">>> Test dft_fftshift");
assert ([0, 1, 2, 3, 4, 5] == dft_fftshift([3, 4, 5, 0, 1, 2]), "Wrong even dft_fftshift");
assert ([0, 1, 2, 3, 4] == dft_fftshift([3, 4, 0, 1, 2]), "Wrong odd dft_fftshift");
assert ([0, 1] == dft_fftshift([1, 0]), "Wrong length two dft_fftshift");
assert ([0] == dft_fftshift([0]), "Wrong length one dft_fftshift");
assert ([] == dft_fftshift([]), "Wrong empty dft_fftshift");


//------------------------------------------------------------------------------
// Test rDFT
//------------------------------------------------------------------------------
N_time = 8;
t_arr = num_linspace(0, 1, N_time);  // N_time points in [0, 1>
K = 2;  // K periods in t_arr for frequency in bin k = K
ampl = 2;  // cos in bin k > 0
dc = 0.1;  // DC in bin k = 0
phi = 30;

// Input signal
x_arr = num_cos_signal(t_arr, K, ampl, phi, dc);

// Expected rDFT(x_arr)
XR_exp = [[0.8, 0, ampl * 3.4641, 0, 0],
          [0, 0, ampl * 2, 0, 0]];

// rDFT
XR = rdft(x_arr);
XR_abs = num_complex_abs(XR[0], XR[1]);
XR_angle = num_complex_angle(XR[0], XR[1]);

echo();
echo(">>> Test rDFT");
echo(t_arr = t_arr);
echo(x_arr = x_arr);
echo("XR_exp real =", XR_exp[0]);
echo("XR_exp imag =", XR_exp[1]);
echo("XR real =", XR[0]);
echo("XR imag =", XR[1]);
echo(XR_abs = XR_abs);
echo(XR_angle = XR_angle);

// . test x_arr cos
assert (num_allclose(XR[0], XR_exp[0], rtol=0, atol=1e-5), "Wrong rDFT real for atol");
assert (num_allclose(XR[1], XR_exp[1], rtol=0, atol=1e-5), "Wrong rDFT imag for atol");
assert (num_allclose(XR[0], XR_exp[0], rtol=1e-6, atol=1e-10), "Wrong rDFT real for rtol");
assert (num_allclose(XR[1], XR_exp[1], rtol=1e-6, atol=1e-10), "Wrong rDFT imag for rtol");

// . test impulse, dc
// . a odd length signal, b even length signal
// . lower case for time domain, upper case for frequency domain
a_impulse_0 = [1, 0, 0, 0, 0, 0, 0];
b_impulse_1 = [0, 1, 0, 0, 0, 0, 0, 0];
a_dc = [1, 1, 1, 1, 1, 1, 1];
b_dc = [1, 1, 1, 1, 1, 1, 1, 1];

A_impulse_0 = rdft(a_impulse_0);
B_impulse_1 = rdft(b_impulse_1);
A_dc = rdft(a_dc);
B_dc = rdft(b_dc);

echo("A_impulse_0 real =", A_impulse_0[0]);
echo("A_impulse_0 imag =", A_impulse_0[1]);
echo("B_impulse_1 real =", B_impulse_1[0]);
echo("B_impulse_1 imag =", B_impulse_1[1]);
echo("A_dc real =", A_dc[0]);
echo("A_dc imag =", A_dc[1]);
echo("B_dc real =", B_dc[0]);
echo("B_dc imag =", B_dc[1]);

assert (num_allclose(A_impulse_0[0], [1, 1, 1, 1]), "Wrong rDFT real A_impulse");
assert (num_allclose(A_impulse_0[1], [0, 0, 0, 0]), "Wrong rDFT imag A_impulse");
assert (num_allclose(B_impulse_1[0], [1, sqrt(0.5), 0, -sqrt(0.5), -1]), "Wrong rDFT real B_impulse");
assert (num_allclose(B_impulse_1[1], [0, -sqrt(0.5), -1, -sqrt(0.5), 0]), "Wrong rDFT imag B_impulse");
assert (num_allclose(A_dc[0], [7, 0, 0, 0]), "Wrong rDFT real A_dc");
assert (num_allclose(A_dc[1], [0, 0, 0, 0]), "Wrong rDFT imag A_dc");
assert (num_allclose(B_dc[0], [8, 0, 0, 0, 0]), "Wrong rDFT real B_dc");
assert (num_allclose(B_dc[1], [0, 0, 0, 0, 0]), "Wrong rDFT imag B_dc");

//------------------------------------------------------------------------------
// All asserts went OK when program gets to here
//------------------------------------------------------------------------------
echo("PASSED");
