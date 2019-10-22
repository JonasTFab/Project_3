# Project_3 Numerical integration

-All the code for this project is contained in the folder "Numerical_integration"

-Test runs of the algoritms can be found in the folder test_run_data

-In folder "Numerical_integration" you will find two scrips.
	"main.cpp" that contains all the algorithms discussed in the report.
	"plots.py" that generates all the plots used in the report 

-The program main.cpp contais various functions which employ different methods of     	numerical integration

	- double integrating_function(double x1, double y1, double z1, double x2, double y2,  double z2)
		-This is the wave function we wish to integrate. Input is coordinates of two particles
		
	-double int_func_spherical_coord(double r1, double r2, double theta1, double theta2, double phi1, double phi2)
		- konverts from kartesian to spherical
	
	-double ran()
		Just a function to make the random initializer look a bit better
	
	-double y(double x){ 
	 	-cumulative function
	-double p(double y){ 
	// radial PDF
	

	-double monte_carlo_improved(int N){
		// The improved Monte Carlo method. This method has the variables changed to
		// spherical coordinate instead of cartesian.
		
	-double brute_monte_carlo(int N, double a, double b){
	// Crude Monte Carlo evaluation with cartesian coordinates to
	// describe the two electrons.
         // N is number of Monte Carlo samples
		
	-void gauss_legendre(double x1, double x2, double x[], double w[], int N)
	
	-void gauss_laguerre(double *x, double *w, int n, double alf)