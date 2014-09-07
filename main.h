#ifndef MAIN_H
#define MAIN_H


double* tridiagonal (double*a,double*b,double*c,double*w, int n);
/* general tridiagonal matrix solver, a is the main diagonal,
 * b and c are the side diagonals, w is the array that contains the results.
 * w is not changed in tridiagonal
 * n is the size of the arrays
 */


double* tridiagonaldiff (double *ah,double* w, int n);
/* specific tridiagonal matrix solver for Poisson equation,
 * results w are not changed in function
 * ah is the static verctor derived from the Gaussian elimination applied to the Poisson matrix
 * n is the size of the arrays
 */


double* f(double x);
/* source function f(x)=100*exp(-10*x)
 */

double* diffreference(double x);

/* reference solution for source function ref(x)=1-(1-exp(-10))*x-exp(-10*x)
 */

void write(double *z, double *y, int n, char *file);
/* writes (z,y) pairs from double vectors of size n into a file
 */

void printmatrix(double ** A, int n, int m);
/* prints a n x m matrix A
 * A is not changed in function
 */

double findmax(double* ei, int n);
/* finds the max value in an array
 * the array is not changed in function
 */

#endif // MAIN_H

