Certified Adaptive-Precision Evaluation of the Voigt Function Using Ball Arithmetic
Overview
Author: Hssaïne AMAMOU
Corresponding Author: Hssaïne AMAMOU <amamou@univ-tln.fr>
Journal: Computer Physics Communications

This repository provides a C++ implementation of a certified adaptive-precision algorithm for evaluating the Voigt function using precision Ball arithmetic.
The software is associated with the scientific manuscript:
Certified adaptive-precision evaluation of the Voigt function using Ball arithmetic
Submitted to Computer Physics Communications.
The implementation uses the FLINT library, including its Arb and Acb precision arithmetic modules.
The program evaluates the Voigt function for user-defined sets of damping parameters and frequency offsets. It automatically increases the working precision until the prescribed relative-accuracy target is achieved or the maximum working precision is reached.
1. Mathematical definition
The Voigt function is evaluated as the real part of the Faddeeva function:
V(aⓜ,u)=Re⁡[w(-u+ia)],

where
w(z)=e^(-z^2 ) erfc⁡(-iz).

Here:
	ais the dimensionless damping parameter;
	uis the dimensionless frequency displacement from the line center;
	w(z)is the Faddeeva function.
The program evaluates
w(-u+ia)=e^(-(-u+ia)^2 ) erfc⁡(a+iu).

The Voigt function is then obtained by extracting the real part.
2. Numerical implementation
The numerical implementation is based on:
	FLINT: Fast Library for Number Theory;
	Arb: precision real Ball arithmetic;
	Acb: precision complex Ball arithmetic.
The evaluation involves the following operations:
	Construction of the complex argument z=-u+ia.
	Evaluation of the exponential term e^(-z^2 ).
	Evaluation of the complex complementary error function erfc⁡(a+iu).
	Multiplication of the two complex quantities.
	Extraction of the real part.
The resulting Arb value provides an enclosure of the computed real quantity.
3. Adaptive-precision strategy
The program uses an automatic precision-increase procedure.
The initial working precision is 256 bits. If the requested relative accuracy is not achieved, the working precision is successively doubled.
The precision sequence is:
Iteration	Working precision
1	256 bits
2	512 bits
3	1024 bits
4	2048 bits
5	4096 bits
6	8192 bits
The calculation is accepted when the Arb relative-accuracy estimate satisfies
"bits obtained"≥180.

If the target accuracy is not reached at 8192 bits, the calculation is reported as unsuccessful.
Important: The reported accuracy is estimated from the Arb enclosure. It is not obtained by comparison with an independently known exact value.
4. Main features
The software provides:
	Evaluation of the Voigt function for multiple values of aand u.
	precision real and complex arithmetic.
	Automatic increase of the working precision.
	A target relative-accuracy estimate of 180 bits.
	A maximum working precision of 8192 bits.
	Direct reading of decimal input strings.
	Display of the calculated Voigt function values.
	Measurement of computation times for each pair (aⓜ,u).
	Display of the total computation time.
The input values are read as decimal strings to avoid unnecessary conversion through standard double-precision floating-point numbers.
5. Requirements
The following software is required:
	A C++ compiler supporting C++17 or later.
	FLINT with Arb and Acb support.
	GMP.
	MPFR.
	CMake (recommended for compilation).
The exact dependency versions should be recorded when reporting computational performance.
6. Compilation
Linux
First, install or build FLINT with the required Arb and Acb functionality.
Then clone the repository:
bash
git clone https://github.com/YOUR_USERNAME/certified-voigt-ball-arithmetic.git

cd certified-voigt-ball-arithmetic
Configure the project:
bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
Compile:
bash
cmake --build build --parallel
The executable is generated in the build directory.
Run the program:
bash
./build/voigt
The exact compilation and linking requirements may depend on the installed FLINT version and system configuration.
7. Input format
The program requests two vectors.
First, enter the number of damping parameters, followed by their values.
For example:
text
5
0
0.5
1
1.5
2
Next, enter the number of frequency offsets, followed by their values.
For example:
text
4
-1
0
1
2
The program evaluates all combinations of the supplied values.
Decimal inputs with many significant digits are accepted, for example:
text
0.123456789012345678901234567890
8. Output
The program displays two matrices.
Matrix 1: Voigt function values
Each row corresponds to a damping parameter a, and each column corresponds to a frequency offset u.
Matrix 2: Computation times
The execution time for each pair (aⓜ,u)is reported in microseconds.
The total computation time is also displayed.
The numerical results are printed with approximately 25 significant decimal digits. This output formatting does not imply that the printed decimal representation contains the full internal precision.
9. Example
An example input file is provided in:
text
examples/example_input.txt
The program can be run interactively by entering the values shown in that file.
10. Reproducibility and performance
The repository is intended to facilitate access to the implementation and reproduction of the numerical calculations presented in the associated manuscript.
For reproducible performance measurements, the following information should be reported:
	Processor model and architecture.
	Operating system.
	C++ compiler and version.
	FLINT version.
	GMP and MPFR versions.
	Compilation optimization options.
	Input values of aand u.
	Working precision reached for each calculation.
	Execution times.
Performance comparisons should be performed under documented and consistent computational conditions.
11. Citation
If you use this software in your research, please cite the associated software release using its Zenodo DOI.
The citation information will be available in the CITATION.cff file and in the Zenodo record.
12. License
This software is distributed under the license specified in the LICENSE file.
13. Contact
For questions concerning the implementation or scientific method, please contact the corresponding author of the associated manuscript.
14. Related publication
Title: Certified adaptive-precision evaluation of the Voigt function using Ball arithmetic.
Journal: Computer Physics Communications.
Status: Submitted.
The publication reference and DOI will be added following publication.

