
/*
===============================================================================
PROGRAM: HIGH-PRECISION EVALUATION OF THE VOIGT FUNCTION V(a,u)
===============================================================================
Author : Hssaïne AMAMOU
1. PROGRAM DESCRIPTION
----------------------

This C++ program calculates the Voigt function V(a,u) for user-defined
sets of damping parameters a and frequency offsets u.

The Voigt function is defined as:

                  V(a,u) = Re[w(-u + i*a)]

where w(z) is the Faddeeva function:

                  w(z) = exp(-z*z) * erfc(-i*z)

and:

                  z = -u + i*a

The Voigt function is widely used in spectroscopy, astrophysics,
plasma physics, radiative transfer, and atomic and molecular spectroscopy.

It describes a spectral line profile resulting from the convolution
of Gaussian and Lorentzian broadening mechanisms.

The two dimensionless variables have the following meanings:

    a : Dimensionless damping parameter.
        It characterizes the relative contribution of Lorentzian
        broadening.

    u : Dimensionless frequency displacement from the line center.
        It represents the spectral position relative to the central
        frequency, expressed in normalized units.

The program evaluates V(a,u) for every combination of the supplied
values of a and u.


2. NUMERICAL METHOD
-------------------

The calculation is performed using the FLINT library and its Arb and
Acb arbitrary-precision arithmetic modules.

The main numerical operations are:

    - Complex multiplication.
    - Complex exponentiation.
    - Complex complementary error function erfc(z).
    - Extraction of the real part of the Faddeeva function.

The Faddeeva function is evaluated using:

                  w(z) = exp(-z*z) * erfc(-i*z)

For the selected complex argument:

                  z = -u + i*a

we obtain:

                  -i*z = a + i*u

Consequently, the program evaluates:

                  w(-u+i*a)
                  = exp(-(-u+i*a)^2) * erfc(a+i*u)

Finally, the Voigt function is obtained from the real part:

                  V(a,u) = Re[w(-u+i*a)]


3. ARBITRARY-PRECISION ARITHMETIC
---------------------------------

Unlike standard double-precision calculations, the Arb library
provides arbitrary-precision real arithmetic based on rigorous
interval arithmetic.

The Acb library provides the corresponding arbitrary-precision
complex arithmetic.

The program uses an automatic precision-increase procedure.

The working precision starts at:

                  256 bits

It is successively doubled:

                  256, 512, 1024, 2048, 4096, 8192 bits

At each precision level, the Voigt function is recalculated.

The relative accuracy of the resulting Arb value is then estimated.

The calculation is accepted when:

                  bits_obtained >= 180

If the target accuracy is not reached at 8192 bits, the calculation
is considered unsuccessful.

IMPORTANT:
The reported relative accuracy is the accuracy estimate provided by
Arb for the resulting interval. It is not a comparison against an
independently known exact value.


4. INPUT DATA
-------------

The program requests two vectors from the user.

First, the number of elements in vector a, followed by the individual
values of a.

Example:

    0
    0.5
    1
    1.5
    2

Second, the number of elements in vector u, followed by the individual
values of u.

Example:

    -1
    0
    1
    2

The input values are stored as strings.

This is intentional: Arb reads the decimal strings directly, avoiding
an unnecessary conversion through the standard double-precision type.

For example, the following decimal input can be supplied:

    0.123456789012345678901234567890


5. CALCULATION PROCEDURE
------------------------

For each combination of a and u, the program:

    1. Reads the two decimal input strings.
    2. Initializes the required Arb and Acb variables.
    3. Sets the complex argument z = -u + i*a.
    4. Calculates z*z.
    5. Calculates exp(-z*z).
    6. Calculates erfc(a+i*u).
    7. Multiplies the two complex quantities.
    8. Extracts the real part of the result.
    9. Estimates the relative accuracy in bits.
   10. Increases the working precision if necessary.

The procedure is repeated until the target accuracy is reached
or the maximum working precision is exceeded.


6. OUTPUT DATA
--------------

The program displays two matrices.

MATRIX 1: VOIGT FUNCTION VALUES

Each row corresponds to one value of a.

Each column corresponds to one value of u.

MATRIX 2: COMPUTATION TIMES

The second matrix contains the execution time, expressed in
microseconds, for each combination of a and u.

The program also displays the total computation time.

The calculation results are printed with 25 significant digits.


7. COMPILATION AND EXECUTION
----------------------------

The program requires:

    - A C++ compiler supporting C++11 or later.
    - The FLINT library.
    - The Arb and Acb modules included in the FLINT installation.

On a Linux system, a typical compilation command is:

    g++ -O2 -std=c++17 voigt.cpp -o voigt -lflint

Depending on the FLINT installation, additional libraries or
linker options may be necessary.

For example:

    g++ -O2 -std=c++17 voigt.cpp -o voigt \
        -lflint -lmpfr -lgmp

Run the program with:

    ./voigt

On Windows, the corresponding FLINT libraries and include directories
must be configured in the compiler or IDE.

The precise linking options depend on the installed FLINT version.


8. IMPORTANT REMARKS
--------------------

- Calculations are performed independently for every pair (a,u).
- The working precision is automatically adjusted.
- The target relative accuracy is 180 bits.
- The maximum working precision is 8192 bits.
- The calculation time depends on the input values and the precision
  required to achieve the target accuracy.
- The input decimal strings are passed directly to Arb.
- All dynamically allocated FLINT objects are released before
  the program terminates.

===============================================================================
END OF PROGRAM DESCRIPTION
===============================================================================
*/


// ============================================================================
// REQUIRED LIBRARIES
// ============================================================================

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <chrono>

#include <flint/flint.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>

using namespace std;
using namespace chrono;


// ============================================================================
// COMPUTATION OF THE VOIGT FUNCTION
// ============================================================================
//
// Calculates:
//
//                  V(a,u) = Re[w(-u + i*a)]
//
// where:
//
//                  w(z) = exp(-z*z) * erfc(-i*z)
//
// The input parameters a and u are provided as decimal strings.
//
// The result is stored in an Arb variable.
//
// The function automatically increases the working precision until
// the required relative accuracy of 180 bits is reached.
//
// Parameters:
//
//   sa                  : Decimal string representing a.
//   su                  : Decimal string representing u.
//   resultat            : Arb variable containing the calculated result.
//   precision_utilisee  : Working precision used for the calculation.
//   bits_obtenus        : Estimated relative accuracy in bits.
//
// Return value:
//
//   true  : The target accuracy has been reached.
//   false : The calculation failed or the target accuracy was not reached.
// ============================================================================

bool calculer_V(const string& sa, const string& su,
                arb_t resultat,
                slong& precision_utilisee,
                slong& bits_obtenus)
{
    arb_t a, u, minus_u;

    acb_t z, z2, minus_z2, minus_iz;
    acb_t exponential, erfc_value, w;

    // Initialize real Arb variables.

    arb_init(a);
    arb_init(u);
    arb_init(minus_u);

    // Initialize complex Acb variables.

    acb_init(z);
    acb_init(z2);
    acb_init(minus_z2);
    acb_init(minus_iz);
    acb_init(exponential);
    acb_init(erfc_value);
    acb_init(w);

    bool succes = false;

    // Automatically increase the working precision.

    for (slong prec = 256; prec <= 8192; prec *= 2)
    {
        // Convert the decimal input strings into Arb numbers.

        if (arb_set_str(a, sa.c_str(), prec) != 0 ||
            arb_set_str(u, su.c_str(), prec) != 0)
        {
            cerr << "Error: invalid numerical input.\n";
            goto nettoyage;
        }

        // Calculate z = -u + i*a.

        arb_neg(minus_u, u);
        acb_set_arb_arb(z, minus_u, a);

        // Calculate the exponential term exp(-z*z).

        acb_mul(z2, z, z, prec);

        acb_neg(minus_z2, z2);

        acb_exp(exponential, minus_z2, prec);

        // Calculate -i*z = a + i*u.

        acb_set_arb_arb(minus_iz, a, u);

        // Calculate the complementary error function erfc(a+i*u).

        acb_hypgeom_erfc(erfc_value, minus_iz, prec);

        // Calculate the Faddeeva function.

        acb_mul(w, exponential, erfc_value, prec);

        // Extract the real part of the Faddeeva function.

        arb_set(resultat, acb_realref(w));

        // Estimate the relative accuracy in bits.

        bits_obtenus = arb_rel_accuracy_bits(resultat);

        // Check whether the required accuracy has been achieved.

        if (bits_obtenus >= 180)
        {
            precision_utilisee = prec;
            succes = true;
            break;
        }
    }

nettoyage:

    // Release all complex Acb variables.

    acb_clear(w);
    acb_clear(erfc_value);
    acb_clear(exponential);
    acb_clear(minus_iz);
    acb_clear(minus_z2);
    acb_clear(z2);
    acb_clear(z);

    // Release all real Arb variables.

    arb_clear(minus_u);
    arb_clear(u);
    arb_clear(a);

    return succes;
}


// ============================================================================
// MAIN PROGRAM
// ============================================================================

int main()
{
    size_t na, nu;

    // =========================================================================
    // 1. INPUT OF THE DAMPING PARAMETER VECTOR a
    // =========================================================================

    cout << "Number of elements in vector a: ";

    if (!(cin >> na) || na == 0)
        return 1;

    vector<string> a(na);

    cout << "Enter the values of a:\n";

    cout << "Example: 0 0.5 1 1.5 2\n";

    for (size_t i = 0; i < na; ++i)
    {
        cout << "a[" << i + 1 << "] = ";

        if (!(cin >> a[i]))
            return 1;
    }


    // =========================================================================
    // 2. INPUT OF THE FREQUENCY OFFSET VECTOR u
    // =========================================================================

    cout << "\nNumber of elements in vector u: ";

    if (!(cin >> nu) || nu == 0)
        return 1;

    vector<string> u(nu);

    cout << "Enter the values of u:\n";

    cout << "Example: -1 0 1 2\n";

    for (size_t j = 0; j < nu; ++j)
    {
        cout << "u[" << j + 1 << "] = ";

        if (!(cin >> u[j]))
            return 1;
    }


    // =========================================================================
    // 3. INITIALIZATION OF RESULT AND COMPUTATION-TIME MATRICES
    // =========================================================================

    // Matrix containing the calculated Voigt function values.

    vector<vector<string>> valeurs(
        na, vector<string>(nu));

    // Matrix containing the computation times in microseconds.

    vector<vector<double>> temps(
        na, vector<double>(nu, 0.0));

    // Initialize the Arb variable used to store each result.

    arb_t resultat;

    arb_init(resultat);

    // Variables used to record the working precision and accuracy.

    slong prec_utilisee = 0;
    slong bits = 0;

    bool tout_ok = true;


    // =========================================================================
    // CONTINUE WITH STEP 2/2
    // =========================================================================   // =========================================================================
    // 4. CALCULATION OF ALL COMBINATIONS (a,u)
    // =========================================================================

    for (size_t i = 0; i < na; ++i)
    {
        for (size_t j = 0; j < nu; ++j)
        {
            // Start the computation timer.

            auto debut = steady_clock::now();

            // Calculate the Voigt function for the current pair (a,u).

            bool succes = calculer_V(
                a[i], u[j], resultat,
                prec_utilisee, bits);

            // Stop the computation timer.

            auto fin = steady_clock::now();

            // Calculate the execution time in microseconds.

            double duree = duration<double, micro>(
                fin - debut).count();

            temps[i][j] = duree;

            // Handle unsuccessful calculations.

            if (!succes)
            {
                tout_ok = false;

                valeurs[i][j] = "FAILURE";

                continue;
            }

            // Convert the Arb result into a decimal string
            // containing approximately 25 significant digits.

            char* v = arb_get_str(resultat, 25, 0);

            valeurs[i][j] = v ? v : "N/A";

            // Release the memory allocated for the output string.

            if (v)
                flint_free(v);
        }
    }


    // =========================================================================
    // 5. DISPLAY OF THE VOIGT FUNCTION VALUE MATRIX
    // =========================================================================

    cout << "\n\n===== VOIGT FUNCTION MATRIX V(a,u) =====\n\n";

    cout << setw(14) << "a \\ u";

    for (size_t j = 0; j < nu; ++j)
        cout << setw(28) << u[j];

    cout << '\n';

    for (size_t i = 0; i < na; ++i)
    {
        cout << setw(14) << a[i];

        for (size_t j = 0; j < nu; ++j)
            cout << setw(28) << valeurs[i][j];

        cout << '\n';
    }


    // =========================================================================
    // 6. DISPLAY OF THE COMPUTATION-TIME MATRIX
    // =========================================================================

    cout << "\n\n===== COMPUTATION TIMES (microseconds) =====\n\n";

    cout << setw(14) << "a \\ u";

    for (size_t j = 0; j < nu; ++j)
        cout << setw(18) << u[j];

    cout << '\n';

    cout << fixed << setprecision(3);

    for (size_t i = 0; i < na; ++i)
    {
        cout << setw(14) << a[i];

        for (size_t j = 0; j < nu; ++j)
            cout << setw(18) << temps[i][j];

        cout << '\n';
    }


    // =========================================================================
    // 7. CALCULATION AND DISPLAY OF THE TOTAL EXECUTION TIME
    // =========================================================================

    double temps_total = 0.0;

    for (size_t i = 0; i < na; ++i)
    {
        for (size_t j = 0; j < nu; ++j)
        {
            temps_total += temps[i][j];
        }
    }

    cout << "\nTotal computation time: "
         << temps_total << " microseconds\n";


    // =========================================================================
    // 8. DISPLAY OF THE PRECISION SETTINGS
    // =========================================================================

    cout << "\nMaximum working precision: 8192 bits.\n";

    cout << "Target relative accuracy: 180 bits.\n";


    // =========================================================================
    // 9. MEMORY CLEANUP
    // =========================================================================

    // Release the Arb result variable.

    arb_clear(resultat);

    // Release FLINT global resources.

    flint_cleanup();


    // =========================================================================
    // 10. PROGRAM TERMINATION
    // =========================================================================

    return tout_ok ? 0 : 1;
}
