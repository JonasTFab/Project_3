#include <iostream>
#include <armadillo>
#include <string>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <fstream>
//using namespace std;
//using namespace arma;
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
} // end of integrating_function


double int_func_spherical_coord(double r1, double r2, double theta1, double theta2, double phi1, double phi2){
    int alpha = 2;
    double cos_beta = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1-phi2);
    double r12 = sqrt(fabs(r1*r1 + r2*r2 - 2*r1*r2*cos_beta));
    if (r12 < 1e-4){
        return 0;
    }
    else{
        double func = exp(-2*alpha*(r1+r2)) / r12;
        return  func;
    }
}   // end of integrating function



double ran(){ //Just a function to make the random initializer look a bit better
    double ran_nr;
    double err = 1e-4;
    double invers_period = RAND_MAX;
    ran_nr = rand()/invers_period;
    if (ran_nr < err){
        return err;
    }
    else if (ran_nr > (1-err)){
        return 1-err;
    }
    else{
        return ran_nr;
    }
} // end of function ran()

double y(double x){ // cumulative function
    return -log(1-x)/4;
}

double p(double y){ // PDF
    return 4*exp(-4*y);
}


// The improved Monte Carlo method. This method has the variables changed to
// spherical coordinate instead of cartesian.
double monte_carlo_improved(int N){
    // N is the number of samples using Monte Carlo
    double x1, x2, F_tilde, variance;
    double r1, r2, theta1, theta2, phi1, phi2;
    double half_jacobi = 4*pow(pi,4);
    double improved_mc = 0;
    double sigma1 = 0;
    double sigma2 = 0;

    // x1 and x2 generates random number between 0 and 1.
    // r, theta and phi is then generated within the respective
    // intervals. Then calculates the integral for each step that will be summarized.
    for (int i=0; i<N; i++){
        x1 = ran();
        x2 = ran();
        r1 = y(x1);
        r2 = y(x2);
        theta1 = pi*ran();
        theta2 = pi*ran();
        phi1 = 2*pi*ran();
        phi2 = 2*pi*ran();
        F_tilde = r1*r1*r2*r2*sin(theta1)*sin(theta2)*int_func_spherical_coord(r1,r2,theta1,theta2,phi1,phi2) / (p(r1)*p(r2));
        improved_mc += F_tilde;
        sigma1 += F_tilde*F_tilde;
    }
    sigma1 /= N;
    sigma2 = (improved_mc/N)*(improved_mc/N);
    variance = sigma1 - sigma2;
    std::cout << "Variance=" << variance << std::endl;
    improved_mc *= half_jacobi / N;
    return improved_mc;
} // end of function mc_improved()



// Crude Monte Carlo evaluation with cartesian coordinates to
// describe the two electrons.
double brute_monte_carlo(int N, double a, double b){
         // N is number of Monte Carlo samples

         double func, variance;
         double x1, y1, z1, x2, y2, z2;
         double crude_mc = 0;
         double sigma1 = 0;
         double sigma2 = 0;
         //arma::vec <double> x_arr = arma::vec(N);

         double jacobi = pow((b-a),6);
         for (int i = 0; i < N; i++){
           //initialize the random numbers
           x1 =  ran()*(b-a)+a;
           y1 =  ran()*(b-a)+a;
           z1 =  ran()*(b-a)+a;
           x2 =  ran()*(b-a)+a;
           y2 =  ran()*(b-a)+a;
           z2 =  ran()*(b-a)+a;
           func = integrating_function(x1,y1,z1,x2,y2,z2);
           crude_mc  += func;
           sigma1 += func*func;
         }
         sigma1 /= N;
         sigma2 = (crude_mc/N)*(crude_mc/N);
         variance = jacobi*(sigma1 - sigma2);
         std::cout << "Variance=" << variance << std::endl;
         crude_mc = crude_mc * jacobi / N;
         return crude_mc;
         } // end of function brute_monte_carlo


// gauss_legendre() + gauss_laguerre() + gammln() = 99% kopi av Morten sin kode

// Function using Gauss-Legendre quadrature to find
// the weight function and the smooth function.
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
} // end of function gauss_legendre()


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


// Gauss-Laguerre function that finds both the weigth
// function and the smooth function.
void gauss_laguerre(double *x, double *w, int n, double alf)
{
    int its;
    double ai,p1,p2,p3,pp,z,z1;

    for (int i=1; i<=n; i++) {
        if (i == 1) {
            z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
        } else if (i == 2) {
            z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
        } else {
            ai=i-2;
            z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
                (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
        }
        for (its=1; its<=max; its++) {
            p1=1.0;
            p2=0.0;
            for (int j=1; j<=n; j++) {
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

void write_to_file(){
  srand(time(NULL));// seed random number generator with the time now
  int num = 9;
  double lamb = 2.3;
  std::ofstream data;
  data.open("brute_monte_carlo_accuracy.txt");
  //data << "N" << "Integral" << "Error" << std::endl;
  for (int i = 1; i < num;i++){
    int N = pow(10,i);
    double mc = brute_monte_carlo(N, -lamb, lamb);
    double error = fabs(mc - 5*pi*pi/(16*16));
    data << N << " " << mc << " " << error << std::endl;
  }
  data.close();
}



int main(){
    int N;
    double lamb;

    //write_to_file();

    std::string method;
    std::cout << "which method (Legendre(le), Laguerre(la), Monte Carlo(mc), improved Monte Carlo(mc_i))? " << std::endl;
    std::cin >> method;
    //method = "mc";


    if (method=="le"){
        std::cout << "Number of integrating points (samples): " << std::endl;
        std::cin >> N;
        std::cout << "Limits of integrations (start = -end. It is found that 2.3 is an ideal value): " << std::endl;
        std::cin >> lamb;

        // variables that is calculated when sending in to the
        // Gauss-Legendre quadrature function.
        double *x = new double [N];
        double *W = new double [N];

        gauss_legendre(-lamb, lamb, x, W, N);
        double int_gauss_legendre = 0;

        // Runs through N^6 operations to estimate the integral
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

        std::cout << "N=" << N << ", lambda=" << lamb << ", I=" << int_gauss_legendre << std::endl << "We want " << 5*pi*pi/(16*16) << std::endl;
        std::cout << "Difference is " << fabs(int_gauss_legendre-5*pi*pi/(16*16)) << std::endl;

    } // end of Legendre method


    else if (method=="la"){
        int alpha = 2;
        std::cout << "Number of integrating points (mesh): " << std::endl;
        std::cin >> N;
        //N = 10;

        // variables that is calculated when sending in to the
        // Gauss-Laguerre quadrature function.
        double *xgl = new double [N+1];
        double *Wgl = new double [N+1];
        double *xt = new double [N];
        double *xp = new double [N];
        double *Wt = new double [N];
        double *Wp = new double [N];
        double alf = 2;
        gauss_legendre(0, pi, xt, Wt, N);
        gauss_legendre(0, 2*pi, xp, Wp, N);
        gauss_laguerre(xgl, Wgl, N, alf);

        // Runs through N^6 operations to estimate the integral
        double int_gauss_laguerre = 0;
        for (int i1 = 1; i1 <= N; i1++){
            for (int i2 = 1; i2 <= N; i2++){
                for (int i3 = 0; i3 < N; i3++){
                    for (int i4 = 0; i4 < N; i4++){
                        for (int i5 = 0; i5 < N; i5++){
                            for (int i6 = 0; i6 < N; i6++){
                                        int_gauss_laguerre += Wgl[i1]*Wgl[i2]*Wt[i3]*Wt[i4]*Wp[i5]*Wp[i6]
                                                *int_func_spherical_coord(xgl[i1],xgl[i2],xt[i3],xt[i4],xp[i5],xp[i6])
                                                *sin(xt[i3])*sin(xt[i4]) / (1024*exp(-2*alpha*(xgl[i1]+xgl[i2])));
                            }
                        }
                    }
                }
            }
        }

        std::cout << "N=" << N << ", I=" << int_gauss_laguerre << std::endl << "We want " << 5*pi*pi/(16*16) << std::endl;
        std::cout << "Difference is " << fabs(int_gauss_laguerre-5*pi*pi/(16*16)) << std::endl;


    } // end of Laguerre method


    else if (method=="mc"){
        std::cout << "Number of integrating points (samples): " << std::endl;
        std::cin >> N;
        std::cout << "Limits of integrations (start = -end): " << std::endl;
        std::cin >> lamb;

        srand(time(NULL));// seed random number generator with the time now
        auto start = std::chrono::high_resolution_clock::now();
        double mc_crude = brute_monte_carlo(N, -lamb, lamb);
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_crude_MC = finish - start;
        std::cout << "N=" << N << ", lambda=" << lamb << ", I=" << mc_crude << std::endl << "We want " << 5*pi*pi/(16*16) << std::endl;
        std::cout << "Difference is " << fabs(mc_crude-5*pi*pi/(16*16)) << std::endl;
        std::cout << "Time brute force Monte carlo: " << elapsed_crude_MC.count() << "s" << std::endl;
    }   // end of brute force Monte Carlo method


    else if (method=="mc_i"){
        std::cout << "Number of integrating points (samples): " << std::endl;
        std::cin >> N;

        srand(time(NULL));// seed random number generator with the time now
        auto start = std::chrono::high_resolution_clock::now();
        double mc_imp = monte_carlo_improved(N);
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_improved_MC = finish - start;
        std::cout << "N=" << N << ", I=" << mc_imp << std::endl << "We want " << 5*pi*pi/(16*16) << std::endl;
        std::cout << "Difference is " << fabs(mc_imp-5*pi*pi/(16*16)) << std::endl;
        std::cout << "Time improved Monte carlo: " << elapsed_improved_MC.count() << "s" << std::endl;



    }   // end of improved Monte Carlo method

    else {
        std::cout << "No valid method! Try again!" << std::endl;
        return 0;
    }



    return 0;
}   // end of main
