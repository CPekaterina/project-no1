#include <iostream>
#include <iomanip>
#include <cmath>
#include "lib.h"
#include <fstream>



using namespace std;

double* tridiagonal (double*a,double*b,double*c,double*w, int n);   //general tridiagonal matrix solver
double* tridiagonaldiff (double *ah,double* w, int n);              //sprecific trdiagonal matrix solver for Poisson equation
double* f(double x);                                                //source function
double* diffreference(double x);                                    //reference solution for source function
void write(double *z, double *y, int n, char *file);                //writes a (z,y) double vector of size n into a file
void printmatrix(double ** A, int n, int m);                        //prints a n x m matrix A


int main()
{
    //get the number of gridpoints

    int n;
    cout << "What is the number of grid points? n= ";
    cin >> n;

    //step length and steparray

    double h =double(1)/double(n+1);
    double* steparray;
    steparray = new double[n];
    for (int i=0;i<n;i++)
    {
        steparray[i]=(i+1)*h;
    }

    //create the right side of the equation

    double*w;
    w = new double[n];
    for (int i=1;i<=n;i++)
    {
        w[i-1]=*f(double(i)*h)*h*h;     // b[i]=f[i]*h, see project description
    }




    //solution with tridiagonal

    double* a,*b,*c;
    a = new double[n];
    b = new double[n];
    c = new double[n];
    c[0]=double(0);
    c[n-1]=double(-1);
    b[0]=double(-1);
    b[n-1]=double(0);
    a[0]=a[n-1]=double(2);
    for (int i=1;i<n-1;i++)
    {
        a[i] = double(2);
        b[i] = double(-1);
        c[i] = double(-1);
    }

    double* results2 = tridiagonal(a,b,c,w,n);


    //solution with tridiagonaldiff

    //w has to be reset because of call by reference if w in tridiagonal
    w = new double[n];
    for (int i=1;i<=n;i++)
    {
        w[i-1]=*f(double(i)*h)*h*h;     // b[i]=f[i]*h, see project description
    }


    //shortcut to save 2*(n-1) flops

    double* ah;
    ah = new double[n];
    ah[0]=2;
    for (int i=1;i<n;i++)
    {
        ah[i] = double(2)-double(1)/ah[i-1];
    }


    //compute the results

    double *results=tridiagonaldiff(ah,w,n);

    //compute the reference solution

    double *ref;
    ref = new double[n];
    for(int i=0;i<n;i++)
    {
        ref[i]=*diffreference(double(1+i)*h);

    }

    //write reference solution in a .dat file and also the results

    char filename[30]={0};
    char reffile[30]={0};
    cout << "Name the results file: ";
    cin >> filename;
    cout << "Name the reference file: ";
    cin >> reffile;

    write(steparray,ref,n,reffile);
    write(steparray,results,n,filename);

    //part D: LU-decomposition

    //write the set of equations as a matrix

    double **A;
    A = new double * [n];
    for (int i = 0; i < n; i++)
    A[i] = new double[n];

    for(int i=0; i<n;i++)
    {
        for(int j=0; j<n;j++)
        {
            A[i][j]=0;
            A[i][i]=double(2);
            A[i][i+1]=double(-1);
            A[i][i-1]=double(-1);
        }
    }
    printmatrix(A, n, n);
    int *indx;
    indx = new int[n];
    double *d;
    ludcmp(A,n,indx,d);

    return 0;
}

// all the functions for LE-solving

double* tridiagonal (double*a,double*b,double*c,double*w, int n)
{
    double* x;
    x = new double[n];

    for (int i=1;i<n;i++)
    {
        double temp = c[i]/a[i-1];
        a[i] -= b[i-1]*temp;
        w[i] -= w[i-1]*temp;
    }

    x[n-1] = w[n-1]/a[n-1];

    for (int i=n-2;i>=0;i--)
    {
        x[i] = (w[i]-b[i]*x[i+1])/a[i];
    }
    return x;
}

double* tridiagonaldiff (double* a, double* w, int n)
{
    double* x;
    x = new double[n];

    for (int i=1;i<n;i++)
    {
        w[i] += w[i-1]/a[i-1];
    }

    x[n-1] = w[n-1]/a[n-1];

    for (int i=n-2;i>=0;i--)
    {
        x[i] = (w[i]+x[i+1])/a[i];
    }
    return x;

}


double* f(double x)
{
    double res = double(100)*exp(double(-10)*x);
    return &res;
}
double* diffreference(double x)
{
    double res = double(1)-(double(1)-exp(-double(10)))*x-exp(double(-10)*x);
    return &res;
}

//function to write results in .dat files

void write(double *z, double *y, int n, char *file)
{
    ofstream resout;
    resout.open(file);
    for (int i=0; i<n; i++)
    {
        resout << setprecision(15) << setw(19) << z[i] << " " << setprecision(15) << setw(19) << y[i] << endl;
    }
    resout.close();
}

void printmatrix(double ** A, int n, int m)
{

    for(int i=0;i<n;i++)
    {
        cout << "| ";
        for(int j=0; j<m; j++)
        {
            cout << setw(6) << A[i][j] << " ";
        }
        cout << "|" << endl;
    }
}
