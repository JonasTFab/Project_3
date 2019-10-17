#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors\n",
           processor_name, world_rank, world_size);

    // Finalize the MPI environment.
    MPI_Finalize();
}

//lager monte carlo greia her for aa slippe aa kjore programmet hele tiden
/*
#include <iostream>
#include <armadillo>
#include <string>
#include "mpi.h"

//new stuff
#include <cstdlib>
#include <ctime>

const double pi = 3.141592653589793238463;
const double eps = 1e-8;
const double eps2 = 3e-14;
const int max = 10;
double exact = 5.*pow(pi,2)/pow(16,2);
using namespace std;
//using namespace arma;
double invers_period = 1./RAND_MAX;

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


int main(int nargs, char* args[]){
  int numprocs, my_rank; // MPI initializations
  MPI_Init (&nargs, &args);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  cout << "Hello world, I have rank " << my_rank << " out of " << numprocs << std::endl; // End MPI
  MPI_Finalize ();
  int n, i;
  n = 1000000; //number of monte carlo samples
  srand(time(NULL));// seed random number generator with the time now
  //Crude monte carlo evaluation
  double crude_mc, sum_sigma, func, variance;
  double x1, y1, z1, x2, y2, z2;
  crude_mc = sum_sigma = 0.;
  for (i = 0; i <= n; i++){

    //initialize the random numbers
    x1 = rand()*invers_period;
    y1 = rand()*invers_period;
    z1 = rand()*invers_period;
    x2 = rand()*invers_period;
    y2 = rand()*invers_period;
    z2 = rand()*invers_period;
    func = integrating_function(x1,y1,z1,x2,y2,z2);
    crude_mc  += func;
    sum_sigma += pow(func,2);
    //std::cout << x1 << x2 << y1 << z1 << std::endl;
  }
  crude_mc = crude_mc/((double) n);
  sum_sigma = sum_sigma/((double) n);
  variance = sum_sigma - crude_mc * crude_mc;

  std::cout << crude_mc<< " " << exact << std::endl;
    return 0;
  }

// g++ monte_carlo.cpp -o mc -larmadillo -llapack -lblas
*/
