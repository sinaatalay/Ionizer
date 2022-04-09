// NOTE BY SINA ATALAY: I have only changed one thing in this code: function mgmres_st won't take
// relative tolerance values as a parameter anymore and won't check it for the two if statements.
// Deleted variable rho_tol.
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "Mgmres.h"
#include "Log.h"

using namespace std;


void atx_st(int n, int nz_num, int ia[], int ja[], double a[], double x[], double w[]) {
    //****************************************************************************80
    //
    //  Purpose:
    //
    //    ATX_ST computes A'*x for a matrix stored in sparse triplet form.
    //
    //  Discussion:
    //
    //    The matrix A is assumed to be sparse.  To save on storage, only
    //    the nonzero entries of A are stored.  For instance, the K-th nonzero
    //    entry in the matrix is stored by:
    //
    //      A(K) = value of entry,
    //      IA(K) = row of entry,
    //      JA(K) = column of entry.
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license. 
    //
    //  Modified:
    //
    //    18 July 2007
    //
    //  Author:
    //
    //    Original C version by Lili Ju.
    //    C++ version by John Burkardt.
    //
    //  Reference:
    //
    //    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
    //    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
    //    Charles Romine, Henk van der Vorst,
    //    Templates for the Solution of Linear Systems:
    //    Building Blocks for Iterative Methods,
    //    SIAM, 1994,
    //    ISBN: 0898714710,
    //    LC: QA297.8.T45.
    //
    //    Tim Kelley,
    //    Iterative Methods for Linear and Nonlinear Equations,
    //    SIAM, 2004,
    //    ISBN: 0898713528,
    //    LC: QA297.8.K45.
    //
    //    Yousef Saad,
    //    Iterative Methods for Sparse Linear Systems,
    //    Second Edition,
    //    SIAM, 2003,
    //    ISBN: 0898715342,
    //    LC: QA188.S17.
    //
    //  Parameters:
    //
    //    Input, int N, the order of the system.
    //
    //    Input, int NZ_NUM, the number of nonzeros.
    //
    //    Input, int IA[NZ_NUM], JA[NZ_NUM], the row and column indices
    //    of the matrix values.
    //
    //    Input, double A[NZ_NUM], the matrix values.
    //
    //    Input, double X[N], the vector to be multiplied by A'.
    //
    //    Output, double W[N], the value of A'*X.
    //

    int i;
    int j;
    int k;

    for (i = 0; i < n; i++)
    {
        w[i] = 0.0;
    }

    for (k = 0; k < nz_num; k++)
    {
        i = ia[k];
        j = ja[k];
        w[j] = w[j] + a[k] * x[i];
    }
    return;
}
//****************************************************************************80

void ax_st(int n, int nz_num, int ia[], int ja[], double a[], double x[], double w[]) {

    //****************************************************************************80
    //
    //  Purpose:
    //
    //    AX_ST computes A*x for a matrix stored in sparse triplet form.
    //
    //  Discussion:
    //
    //    The matrix A is assumed to be sparse.  To save on storage, only
    //    the nonzero entries of A are stored.  For instance, the K-th nonzero
    //    entry in the matrix is stored by:
    //
    //      A(K) = value of entry,
    //      IA(K) = row of entry,
    //      JA(K) = column of entry.
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license. 
    //
    //  Modified:
    //
    //    18 July 2007
    //
    //  Author:
    //
    //    Original C version by Lili Ju.
    //    C++ version by John Burkardt.
    //
    //  Reference:
    //
    //    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
    //    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
    //    Charles Romine, Henk van der Vorst,
    //    Templates for the Solution of Linear Systems:
    //    Building Blocks for Iterative Methods,
    //    SIAM, 1994,
    //    ISBN: 0898714710,
    //    LC: QA297.8.T45.
    //
    //    Tim Kelley,
    //    Iterative Methods for Linear and Nonlinear Equations,
    //    SIAM, 2004,
    //    ISBN: 0898713528,
    //    LC: QA297.8.K45.
    //
    //    Yousef Saad,
    //    Iterative Methods for Sparse Linear Systems,
    //    Second Edition,
    //    SIAM, 2003,
    //    ISBN: 0898715342,
    //    LC: QA188.S17.
    //
    //  Parameters:
    //
    //    Input, int N, the order of the system.
    //
    //    Input, int NZ_NUM, the number of nonzeros.
    //
    //    Input, int IA[NZ_NUM], JA[NZ_NUM], the row and column indices
    //    of the matrix values.
    //
    //    Input, double A[NZ_NUM], the matrix values.
    //
    //    Input, double X[N], the vector to be multiplied by A.
    //
    //    Output, double W[N], the value of A*X.
    //
    int i;
    int j;
    int k;

    for (i = 0; i < n; i++)
    {
        w[i] = 0.0;
    }

    for (k = 0; k < nz_num; k++)
    {
        i = ia[k];
        j = ja[k];
        w[i] = w[i] + a[k] * x[j];
    }
    return;
}
//****************************************************************************80

void mgmres_st(int n, int nz_num, int ia[], int ja[], double a[], double x[], double rhs[], int itr_max, int mr, double tol_abs)

//****************************************************************************80
//
//  Purpose:
//
//    MGMRES_ST applies restarted GMRES to a matrix in sparse triplet form.
//
//  Discussion:
//
//    The linear system A*X=B is solved iteratively.
//
//    The matrix A is assumed to be stored in sparse triplet form.  Only
//    the nonzero entries of A are stored.  For instance, the K-th nonzero
//    entry in the matrix is stored by:
//
//      A(K) = value of entry,
//      IA(K) = row of entry,
//      JA(K) = column of entry.
//
//    The "matrices" H and V are treated as one-dimensional vectors
//    which store the matrix data in row major form.
//
//    This requires that references to H[I][J] be replaced by references
//    to H[I+J*(MR+1)] and references to V[I][J] by V[I+J*N].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 July 2007
//
//  Author:
//
//    Original C version by Lili Ju.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
//    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
//    Charles Romine, Henk van der Vorst,
//    Templates for the Solution of Linear Systems:
//    Building Blocks for Iterative Methods,
//    SIAM, 1994,
//    ISBN: 0898714710,
//    LC: QA297.8.T45.
//
//    Tim Kelley,
//    Iterative Methods for Linear and Nonlinear Equations,
//    SIAM, 2004,
//    ISBN: 0898713528,
//    LC: QA297.8.K45.
//
//    Yousef Saad,
//    Iterative Methods for Sparse Linear Systems,
//    Second Edition,
//    SIAM, 2003,
//    ISBN: 0898715342,
//    LC: QA188.S17.
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Input, int NZ_NUM, the number of nonzero matrix values.
//
//    Input, int IA[NZ_NUM], JA[NZ_NUM], the row and column indices
//    of the matrix values.
//
//    Input, double A[NZ_NUM], the matrix values.
//
//    Input/output, double X[N]; on input, an approximation to
//    the solution.  On output, an improved approximation.
//
//    Input, double RHS[N], the right hand side of the linear system.
//
//    Input, int ITR_MAX, the maximum number of (outer) iterations to take.
//
//    Input, int MR, the maximum number of (inner) iterations to take.
//    MR must be less than N.
//
//    Input, double TOL_ABS, an absolute tolerance applied to the
//    current residual.
//
//    Input, double TOL_REL, a relative tolerance comparing the
//    current residual to the initial residual.
//
{
    double av;
    double* c;
    double delta = 1.0e-03;
    double* g;
    double* h;
    double htmp;
    int i;
    int itr;
    int itr_used;
    int j;
    int k;
    int k_copy;
    double mu;
    double* r;
    double rho;
    double* s;
    double* v;
    bool verbose = true;
    double* y;

    c = new double[mr];
    g = new double[mr + 1];
    h = new double[(mr + 1) * mr];
    r = new double[n];
    s = new double[mr];
    v = new double[n * (mr + 1)];
    y = new double[mr + 1];

    itr_used = 0;

    if (n < mr)
    {
        cerr << "\n";
        cerr << "MGMRES_ST - Fatal error!\n";
        cerr << "  N < MR.\n";
        cerr << "  N = " << n << "\n";
        cerr << "  MR = " << mr << "\n";
        exit(1);
    }

    for (itr = 1; itr <= itr_max; itr++)
    {
        ax_st(n, nz_num, ia, ja, a, x, r);

        for (i = 0; i < n; i++)
        {
            r[i] = rhs[i] - r[i];
        }

        rho = sqrt(r8vec_dot(n, r, r));

        if (verbose)
        {
            cout << "  ITR = " << itr << "  Residual = " << rho << "\n";
        }

        for (i = 0; i < n; i++)
        {
            v[i + 0 * n] = r[i] / rho;
        }

        g[0] = rho;
        for (i = 1; i <= mr; i++)
        {
            g[i] = 0.0;
        }

        for (i = 0; i < mr + 1; i++)
        {
            for (j = 0; j < mr; j++)
            {
                h[i + j * (mr + 1)] = 0.0;
            }
        }

        for (k = 1; k <= mr; k++)
        {
            k_copy = k;

            ax_st(n, nz_num, ia, ja, a, v + (k - 1) * n, v + k * n);

            av = sqrt(r8vec_dot(n, v + k * n, v + k * n));

            for (j = 1; j <= k; j++)
            {
                h[(j - 1) + (k - 1) * (mr + 1)] = r8vec_dot(n, v + k * n, v + (j - 1) * n);
                for (i = 0; i < n; i++)
                {
                    v[i + k * n] = v[i + k * n] - h[(j - 1) + (k - 1) * (mr + 1)] * v[i + (j - 1) * n];
                }
            }

            h[k + (k - 1) * (mr + 1)] = sqrt(r8vec_dot(n, v + k * n, v + k * n));

            if ((av + delta * h[k + (k - 1) * (mr + 1)]) == av)
            {
                for (j = 1; j <= k; j++)
                {
                    htmp = r8vec_dot(n, v + k * n, v + (j - 1) * n);
                    h[(j - 1) + (k - 1) * (mr + 1)] = h[(j - 1) + (k - 1) * (mr + 1)] + htmp;
                    for (i = 0; i < n; i++)
                    {
                        v[i + k * n] = v[i + k * n] - htmp * v[i + (j - 1) * n];
                    }
                }
                h[k + (k - 1) * (mr + 1)] = sqrt(r8vec_dot(n, v + k * n, v + k * n));
            }

            if (h[k + (k - 1) * (mr + 1)] != 0.0)
            {
                for (i = 0; i < n; i++)
                {
                    v[i + k * n] = v[i + k * n] / h[k + (k - 1) * (mr + 1)];
                }
            }

            if (1 < k)
            {
                for (i = 1; i <= k + 1; i++)
                {
                    y[i - 1] = h[(i - 1) + (k - 1) * (mr + 1)];
                }
                for (j = 1; j <= k - 1; j++)
                {
                    mult_givens(c[j - 1], s[j - 1], j - 1, y);
                }
                for (i = 1; i <= k + 1; i++)
                {
                    h[i - 1 + (k - 1) * (mr + 1)] = y[i - 1];
                }
            }
            mu = sqrt(pow(h[(k - 1) + (k - 1) * (mr + 1)], 2)
                + pow(h[k + (k - 1) * (mr + 1)], 2));
            c[k - 1] = h[(k - 1) + (k - 1) * (mr + 1)] / mu;
            s[k - 1] = -h[k + (k - 1) * (mr + 1)] / mu;
            h[(k - 1) + (k - 1) * (mr + 1)] = c[k - 1] * h[(k - 1) + (k - 1) * (mr + 1)]
                - s[k - 1] * h[k + (k - 1) * (mr + 1)];
            h[k + (k - 1) * (mr + 1)] = 0;
            mult_givens(c[k - 1], s[k - 1], k - 1, g);

            rho = fabs(g[k]);

            itr_used = itr_used + 1;

            if (verbose)
            {
                cout << "  K =   " << k << "  Residual = " << rho << "\n";
            }

            if (rho <= tol_abs)
            {
                break;
            }
        }

        k = k_copy - 1;
        y[k] = g[k] / h[k + k * (mr + 1)];

        for (i = k; 1 <= i; i--)
        {
            y[i - 1] = g[i - 1];
            for (j = i + 1; j <= k + 1; j++)
            {
                y[i - 1] = y[i - 1] - h[(i - 1) + (j - 1) * (mr + 1)] * y[j - 1];
            }
            y[i - 1] = y[i - 1] / h[(i - 1) + (i - 1) * (mr + 1)];
        }

        for (i = 1; i <= n; i++)
        {
            for (j = 1; j <= k + 1; j++)
            {
                x[i - 1] = x[i - 1] + v[(i - 1) + (j - 1) * n] * y[j - 1];
            }
        }

        if (rho <= tol_abs)
        {
            break;
        }
    }

    if (verbose)
    {
        cout << "\n";
        cout << "MGMRES_ST\n";
        cout << "  Number of iterations = " << itr_used << "\n";
        cout << "  Final residual = " << rho << "\n";
    }

    //
    //  Free memory.
    //
    delete[] c;
    delete[] g;
    delete[] h;
    delete[] r;
    delete[] s;
    delete[] v;
    delete[] y;

    return;
}
//****************************************************************************80

void mult_givens(double c, double s, int k, double g[])

//****************************************************************************80
//
//  Purpose:
//
//    MULT_GIVENS applies a Givens rotation to two successive entries of a vector.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 August 2006
//
//  Author:
//
//    Original C version by Lili Ju.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
//    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
//    Charles Romine, Henk van der Vorst,
//    Templates for the Solution of Linear Systems:
//    Building Blocks for Iterative Methods,
//    SIAM, 1994,
//    ISBN: 0898714710,
//    LC: QA297.8.T45.
//
//    Tim Kelley,
//    Iterative Methods for Linear and Nonlinear Equations,
//    SIAM, 2004,
//    ISBN: 0898713528,
//    LC: QA297.8.K45.
//
//    Yousef Saad,
//    Iterative Methods for Sparse Linear Systems,
//    Second Edition,
//    SIAM, 2003,
//    ISBN: 0898715342,
//    LC: QA188.S17.
//
//  Parameters:
//
//    Input, double C, S, the cosine and sine of a Givens
//    rotation.
//
//    Input, int K, indicates the location of the first vector entry.
//
//    Input/output, double G[K+2], the vector to be modified.  On output,
//    the Givens rotation has been applied to entries G(K) and G(K+1).
//
{
    double g1;
    double g2;

    g1 = c * g[k] - s * g[k + 1];
    g2 = s * g[k] + c * g[k + 1];

    g[k] = g1;
    g[k + 1] = g2;

    return;
}
//****************************************************************************80

double r8vec_dot(int n, double a1[], double a2[])

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT computes the dot product of a pair of R8VEC's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT, the dot product of the vectors.
//
{
    int i;
    double value;

    value = 0.0;
    for (i = 0; i < n; i++)
    {
        value = value + a1[i] * a2[i];
    }
    return value;
}
//****************************************************************************80

double* r8vec_uniform_01(int n, int& seed)

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2^31 - 1 )
//      unif = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01[N], the vector of pseudorandom values.
//
{
    int i;
    int k;
    double* r;

    if (seed == 0)
    {
        cerr << "\n";
        cerr << "R8VEC_UNIFORM_01 - Fatal error!\n";
        cerr << "  Input value of SEED = 0.\n";
        exit(1);
    }

    r = new double[n];

    for (i = 0; i < n; i++)
    {
        k = seed / 127773;

        seed = 16807 * (seed - k * 127773) - k * 2836;

        if (seed < 0)
        {
            seed = seed + 2147483647;
        }

        r[i] = (double)(seed) * 4.656612875E-10;
    }

    return r;
}
//****************************************************************************80