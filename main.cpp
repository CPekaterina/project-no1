#include <iostream>
#include <iomanip>
#include <cmath>
#include <lib.h>
#include <fstream>


using namespace std;

double* tridiagonal (double*a,double*b,double*c,double*w, int n);   //general tridiagonal matrix solver
double* tridiagonaldiff (double *ah,double* w, int n);              //specific tridiagonal matrix solver for Poisson equation
double* f(double x);                                                //source function
double* diffreference(double x);                                    //reference solution for source function
<<<<<<< HEAD
double findmax(double* ei, int n);                                  //finds the max value in an array
=======
void write(double *z, double *y, int n, char *file);

>>>>>>> a7cc8d6c107174a87bd18168845af3f4853023cb

int main()
{

    //get the number of gridpoints

    int n;
    cout << "What is the number of grid points? n= ";
    cin >> n;

<<<<<<< HEAD

    //step length
=======
    //step length and steparray
>>>>>>> a7cc8d6c107174a87bd18168845af3f4853023cb

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
        w[i-1]=*f(double(i)*h)*h*h;     // b[i]=f[i]*h*h, see project description
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
        w[i-1]=*f(double(i)*h)*h*h;     // b[i]=f[i]*h*h, see project description
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

<<<<<<< HEAD

    //comparation with reference solution
=======
    //compute the reference solution
>>>>>>> a7cc8d6c107174a87bd18168845af3f4853023cb

    double *ref;
    ref = new double[n];
    for(int i=0;i<n;i++)
    {
        ref[i]=*diffreference(double(1+i)*h);

    }

<<<<<<< HEAD

    //calculate the max value of the relative error

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
    cout << logh << " " << maxvalue << " " << n << endl;
=======
    //write reference solution in a .dat file and also the results

    char filename[30]={0};
    char reffile[30]={0};
    cout << "Name the results file: ";
    cin >> filename;
    cout << "Name the reference file: ";
    cin >> reffile;

    write(steparray,ref,n,reffile);
    write(steparray,results,n,filename);
>>>>>>> a7cc8d6c107174a87bd18168845af3f4853023cb

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

<<<<<<< HEAD

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
=======
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
>>>>>>> a7cc8d6c107174a87bd18168845af3f4853023cb
}
