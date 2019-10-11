#include <iostream>
#include <armadillo>

//using namespace std;
//using namespace arma;

double integrating_function(arma::vec r1_vec, arma::vec r2_vec){
    int alpha, n;
    double func, r1, r2, r_diff;
    arma::Col<double> r_diff_vec;

    alpha = 2;
    n = r1_vec.n_rows;
    r1=0; r2=0; r_diff=0;
    r_diff_vec = r1_vec - r2_vec;
    for(int i=0; i<n; i++){
        r1 += pow(r1_vec(i),2);
        r2 += pow(r2_vec(i),2);
        r_diff += pow(r_diff_vec(i),2);
    }
    r1 = sqrt(r1); r2 = sqrt(r2); r_diff = sqrt(r_diff);

    func = exp(-2*alpha*(r1+r2)) / r_diff;
    return func;
}   // end of integrating function


double legendre_polynomial(int k, double x){
    double l_plus,l_minus,l;
    l = 0;
    l_plus = 1;
    // recursion relation
    for(int i=0; i<k; i++){
        l_minus = l;
        l = l_plus;
        l_plus = (-i*l_minus + (2*i+1)*x*l) / (i+1);
    }
    return l_plus;
}   // end of Legendre polynomial function



//void Gauss_Lagrange(int N){
//    int pol_deg;
//    pol_deg = 2*N-1;
//    roots
//}



int main(){
    int k;
    double x;
    k = 2;
    x = 2;
    legendre_polynomial(k,x);

    return 0;
}


//int main(){
//    double lambda;
//    int n;
//    std::cout << "Number of integrating points (mesh): " << std::endl;
//    std::cin >> n;
//    std::cout << "Limits of integrations (start = -end): " << std::endl;
//    std::cin >> lambda;

//    double x, W;

//    gauss_legendre(-lambda,lambda,x,W,n);

//    return 0;
//}   // end of main
