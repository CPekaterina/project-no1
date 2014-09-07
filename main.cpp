#include <iostream>
#include <iomanip>
#include <cmath>
#include "windows.h"
#include "lib.h"
#include <fstream>
#include "time.h"
#include "main.h"


using namespace std;




int main()
{

    cout << "Solvers for the Poisson differential equation" << endl;
    cout << "---------------------------------------------" << endl;

 //get the number of gridpoints

    int n;
    cout << "What is the number of grid points? n= ";
    cin >> n;
    cout << endl;

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
        w[i-1]=*f(double(i)*h)*h*h; // b[i]=f[i]*h*h, see project description
    }




 //solution with tridiagonal

 //set up the diagonals

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

 //shortcut to save 2*(n-1) flops

    double* ah;
    ah = new double[n];
    ah[0]=2;
    for (int i=1;i<n;i++)
    {
        ah[i] = double(2)-double(1)/ah[i-1];
    }




    cout << endl << "THE TRIDIAGONAL SOLVER (tailored to the Poisson equation)" << endl << endl;



 //compute the results and measure the computation time


    clock_t start, finish;// declare start and final time
    start = clock();

    double *results=tridiagonaldiff(ah,w,n);
    finish = clock();
    clock_t t=finish-start;

    cout << "Execution time for results with tridiagonaldiff: " <<setprecision(10) << ((double)t)/CLOCKS_PER_SEC << " sec." << endl;





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
    cout << endl;

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



 //use ludcmp and lubksb to solve the matrix equation with the right side w[]


    cout << endl << "THE LUDCMP AND LUBKSB SOLVERS (lib.cpp)" << endl << endl;

    start = clock();

    int *indx;
    indx = new int[n];
    double *d;
    ludcmp(A,n,indx,d);
    lubksb(A,n,indx,w);

    finish= clock();
    t=finish-start;
    cout << "Execution time for results with ludcmp and lubksb: "<< setprecision(10) << ((double)t)/CLOCKS_PER_SEC << " sec." << endl;



 //write reference solution in a .dat file and also the results

    char LUname[30]={0};
    cout << "Name the results file for LU: ";
    cin >> LUname;

    write(steparray,results,n,LUname);



 //part C

 //calculate the max value of the relative error

    cout << endl << endl << "CALCULATION OF THE MAXIMUM RELATIVE ERROR" << endl << endl;

    double* ei;
    ei = new double[n];

    for(int i=0;i<n;i++)
    {
      double ref = *diffreference(double(1+i)*h);
      ei[i] = log10(fabs((results[i]- ref)/ref));
    }

    double maxvalue = findmax(ei, n);
    double logh = log10(h);

    cout.precision(10);
    cout << "The max relative error is to find at log(h)=" << logh << " and its value is " << maxvalue << " for n=" << n << "." << endl << endl;

    return 0;
}

// all the functions for LE-solving

double* tridiagonal (double*a,double*b,double*c,double*w, int n)
{
    double* x;
    x = new double[n];
    double* tw;
    tw = new double[n];
    tw[0]=w[0];

    for (int i=1;i<n;i++)
    {
        double temp = c[i]/a[i-1];
        a[i] -= b[i-1]*temp;
        tw[i] = w[i]-tw[i-1]*temp;
    }

    x[n-1] = tw[n-1]/a[n-1];

    for (int i=n-2;i>=0;i--)
    {
        x[i] = (tw[i]-b[i]*x[i+1])/a[i];
    }
    return x;
}

double* tridiagonaldiff (double* a, double* w, int n)
{
    double* x;
    x = new double[n];
    double* tw;
    tw = new double[n];
    tw[0]=w[0];

    for (int i=1;i<n;i++)
    {
        tw[i] =w[i]+tw[i-1]/a[i-1];
    }

    x[n-1] = tw[n-1]/a[n-1];

    for (int i=n-2;i>=0;i--)
    {
        x[i] = (tw[i]+x[i+1])/a[i];
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

//function to find the maximum of an array ei of size n

double findmax(double* ei, int n)
{
    double max = ei[0];

    for(int i=1;i<n; i++)
    {
        if(ei[i]>max)
        {
            max = ei[i];
        }
    }
    return max;
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

//function to print a matrix as a check for LU

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
