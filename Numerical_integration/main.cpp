#include <iostream>
#include <armadillo>

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



// gauleg = 99.9% kopi av Morten sin kode
       /*
       ** The function
       **              gauleg()
       ** takes the lower and upper limits of integration x1, x2, calculates
       ** and return the abcissas in x[0,...,n - 1] and the weights in w[0,...,n - 1]
       ** of length n of the Gauss--Legendre n--point quadrature formulae.
       */

void gauleg(double x1, double x2, double x[], double w[], int N)
{
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1,eps;
   double      const  pi = 3.14159265359;
   double      *x_low, *x_high, *w_low, *w_high;
   eps = 1e-8;

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
} // End_ function gauleg()



int main(){
    double lamb;
    int N;
    //std::cout << "Number of integrating points (mesh): " << std::endl;
    //std::cin >> N;
    //std::cout << "Limits of integrations (start = -end): " << std::endl;
    //std::cin >> lamb;
    lamb = 2.3;     // given from the python plot of the electron
    N = 25;
    const double pi = 3.141592653589793238463;

    double *x = new double [N];
    double *W = new double [N];
    gauleg(-lamb, lamb, x, W, N);


    double int_gauss = 0;
    for (int i1 = 0; i1 < N; i1++){
        for (int i2 = 0; i2 < N; i2++){
            for (int i3 = 0; i3 < N; i3++){
                for (int i4 = 0; i4 < N; i4++){
                    for (int i5 = 0; i5 < N; i5++){
                        for (int i6 = 0; i6 < N; i6++){
                                    int_gauss += W[i1]*W[i2]*W[i3]*W[i4]*W[i5]*W[i6]*
                                            integrating_function(x[i1],x[i2],x[i3],x[i4],x[i5],x[i6]);
                        }
                    }
                }
            }
        }
    }

    std::cout << int_gauss << std::endl << "We want " << 5*pi*pi/(16*16) << std::endl;
    std::cout << "Difference is " << fabs(int_gauss-5*pi*pi/(16*16)) << std::endl;


    return 0;
}   // end of main
