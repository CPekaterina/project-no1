#include <iostream>
#include<cmath>

using namespace std;

double* tridiagonal (double*a,double*b,double*c,double*w, int n);   //general tridiagonal matrix solver
double* tridiagonaldiff (double* w, int n);                         //sprecific trdiagonal matrix solver for Poisson equation
double* f(double x);                                                //source function
double* diffreference(double x);                                    //reference solution for source function

int main()
{
    //get the number of gridpoints

    int n;
    cout << "What is the number of grid points? n= ";
    cin >> n;

    //step length

    double h =double(1)/double(n+1);

    //create the right side of the equation

    double*w;
    w = new double[n];
    for (int i=1;i<=n;i++)
    {
        w[i-1]=*f(double(i)*h)*h*h;     // b[i]=f[i]*h, see project description
    }

    //compute the results

    double *results=tridiagonaldiff(w,n);

    //comparation with reference solution

    for(int i=0;i<n;i++)
    {
        cout << results[i] << " " << *diffreference(double(1+i)*h) << endl;
    }

    return 0;
}

// all the functions for LE-solving

double* tridiagonal (double*a,double*b,double*c,double*w, int n)
{
    double* x, *ah, *wh;
    x = new double[n];
    ah = new double[n];
    wh = new double[n];

    ah[0] =a[0];
    wh[0] =w[0];

    for (int i=1;i<n;i++)
    {
        double temp = c[i]/a[i-1];
        ah[i] = a[i]-b[i-1]*temp;
        wh[i] = w[i]-wh[i-1]*temp;
    }

    x[n-1] = wh[n-1]/ah[n-1];
    for (int i=n-2;i>=0;i--)
    {
        x[i] = (wh[i]-b[i+1]*x[i+1])/ah[i];
    }
    return x;
}


double* tridiagonaldiff (double* w, int n)
{
    double* x, *ah, *wh;
    x = new double[n];
    ah = new double[n];
    wh = new double[n];

    ah[0] = 2;
    wh[0] = w[0];

    for (int i=1;i<n;i++)
    {
        ah[i] = double(2)-double(1)/ah[i-1];
        wh[i] = w[i]+wh[i-1]/ah[i-1];
    }

    x[n-1] = wh[n-1]/ah[n-1];

    for (int i=n-2;i>=0;i--)
    {
        x[i] = (wh[i]+x[i+1])/ah[i];
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
