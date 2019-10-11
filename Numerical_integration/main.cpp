#include <iostream>
#include <armadillo>
#include <string>

const double pi = 3.141592653589793238463;
const double eps = 1e-8;
const double eps2 = 3e-14;
const int max = 10;


//using namespace std;
//using namespace arma;

double integrating_function(double x1, double y1, double z1, double x2, double y2, double z2){
    int alpha = 2;
    double r1, r2, r_diff;

    r1 = sqrt(x1*x1+y1*y1+z1*z1);
    r2 = sqrt(x2*x2+y2*y2+z2*z2);
    r_diff = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
    if (r_diff < 1e-4){
        return 0;
    }
    else{
        return exp(-2*alpha*(r1+r2)) / r_diff;
    }
}   // end of integrating function



// gauss_legendre() + gauss_laguerre() +  gammln() = 99% kopi av Morten sin kode
       /*
       ** The function
       **              gauleg()
       ** takes the lower and upper limits of integration x1, x2, calculates
       ** and return the abcissas in x[0,...,n - 1] and the weights in w[0,...,n - 1]
       ** of length n of the Gauss--Legendre n--point quadrature formulae.
       */

void gauss_legendre(double x1, double x2, double x[], double w[], int N)
{
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      *x_low, *x_high, *w_low, *w_high;

   m  = (N + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;                                       // pointer initialization
   x_high = x + N - 1;
   w_low  = w;
   w_high = w + N - 1;

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(N + 0.5));

           /*
           ** Starting with the above approximation to the ith root
           ** we enter the mani loop of refinement bt Newtons method.
           */

      do {
        p1 = 1.0;
        p2 = 0.0;

       /*
       ** loop up recurrence relation to get the
       ** Legendre polynomial evaluated at x
       */

     for(j = 1; j <= N; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
     }

       /*
       ** p1 is now the desired Legrendre polynomial. Next compute
       ** ppp its derivative by standard relation involving also p2,
       ** polynomial of one lower order.
       */

     pp = N * (z * p1 - p2)/(z * z - 1.0);
     z1 = z;
     z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > eps);

          /*
          ** Scale the root to the desired interval and put in its symmetric
          ** counterpart. Compute the weight and its symmetric counterpart
          */

      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
   }
} // End of function gauss_legendre()


double gammln(double xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
} // end of function gammln()


void gauss_laguerre(double *x, double *w, int n, double alf)
{
    int i,its,j;
    double ai;
    double p1,p2,p3,pp,z,z1;


    for (i=1;i<=n;i++) {
        if (i == 1) {
            z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
        } else if (i == 2) {
            z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
        } else {
            ai=i-2;
            z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
                (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
        }
        for (its=1;its<=max;its++) {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=n;j++) {
                p3=p2;
                p2=p1;
                p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
            }
            pp=(n*p1-(n+alf)*p2)/z;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= eps2) break;
        }
        if (its > max) std::cout << "too many iterations in gauss_laguerre" << std::endl;
        x[i]=z;
        w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
    }
} // end of function gauss_legendre()




int main(){
    double lamb;
    int N;
    std::string method;

    //std::cout << "which method (Legendre(le), Laguerre(la)? " << std::endl;
    //std::cin >> method;
    method = "la";


    if (method=="le"){
        //std::cout << "Number of integrating points (mesh): " << std::endl;
        //std::cin >> N;
        //std::cout << "Limits of integrations (start = -end): " << std::endl;
        //std::cin >> lamb;
        lamb = 2.3;     // given from the python plot of the electron
        N = 25;
        double *x = new double [N];
        double *W = new double [N];
        gauss_legendre(-lamb, lamb, x, W, N);

        double int_gauss_legendre = 0;
        for (int i1 = 0; i1 < N; i1++){
            for (int i2 = 0; i2 < N; i2++){
                for (int i3 = 0; i3 < N; i3++){
                    for (int i4 = 0; i4 < N; i4++){
                        for (int i5 = 0; i5 < N; i5++){
                            for (int i6 = 0; i6 < N; i6++){
                                        int_gauss_legendre += W[i1]*W[i2]*W[i3]*W[i4]*W[i5]*W[i6]*
                                                integrating_function(x[i1],x[i2],x[i3],x[i4],x[i5],x[i6]);
                            }
                        }
                    }
                }
            }
        }

        std::cout << int_gauss_legendre << std::endl << "We want " << 5*pi*pi/(16*16) << std::endl;
        std::cout << "Difference is " << fabs(int_gauss_legendre-5*pi*pi/(16*16)) << std::endl;

    } // end of Legendre method


    else if (method=="la"){
        //std::cout << "Number of integrating points (mesh): " << std::endl;
        //std::cin >> N;
        N = 25;

        double *xgl = new double [N+1];
        double *Wgl = new double [N+1];
        double alf = 1;
        gauss_laguerre(xgl, Wgl, N, alf);

        double int_gauss_laguerre = 0;
        for (int i1 = 1; i1 <= N; i1++){
            for (int i2 = 1; i2 <= N; i2++){
                for (int i3 = 1; i3 <= N; i3++){
                    for (int i4 = 1; i4 <= N; i4++){
                        for (int i5 = 1; i5 <= N; i5++){
                            for (int i6 = 1; i6 <= N; i6++){
                                        int_gauss_laguerre += Wgl[i1]*Wgl[i2]*Wgl[i3]*Wgl[i4]*Wgl[i5]*Wgl[i6]*
                                                integrating_function(xgl[i1],xgl[i2],xgl[i3],xgl[i4],xgl[i5],xgl[i6]);
                            }
                        }
                    }
                }
            }
        }

        std::cout << int_gauss_laguerre << std::endl << "We want " << 5*pi*pi/(16*16) << std::endl;
        std::cout << "Difference is " << fabs(int_gauss_laguerre-5*pi*pi/(16*16)) << std::endl;


    } // end of Laguerre method


    else {
        std::cout << "No valid method! Try again!" << std::endl;
        return 0;
    }








    return 0;
}   // end of main
