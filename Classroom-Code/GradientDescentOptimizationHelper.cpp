#include <Rcpp.h>
#include <cstdlib>
#include <cmath>
#include <list>
#include<math.h> 
using namespace Rcpp;

double e(double x, double y, int k = 120){
  return (1 - x/k)*(1 - y/k)*(1 + x/k)*(1 + y/k)*(-(y + 47) * sin(sqrt(abs(y + x/2 + 47))) - x*sin(sqrt(abs(x - (y + 47)))));
}

std::vector<double> e_prime(std::vector<double> par, int k = 120){
  double x = par[0];
  double y = par[1];
  double dy = -((1-x/k)*(x/k+1)*(y/k+1)*((-y-47)*sin(sqrt(abs(y+x/2+47)))-x*sin(sqrt(abs(y-x+47)))))/k+((1-x/k)*(x/k+1)*(1-y/k)*((-y-47)*sin(sqrt(abs(y+x/2+47)))-x*sin(sqrt(abs(y-x+47)))))/k+(1-x/k)*(x/k+1)*(1-y/k)*(y/k+1)*(-sin(sqrt(abs(y+x/2+47)))+((-y-47)*(y+x/2+47)*cos(sqrt(abs(y+x/2+47))))/(2*pow(abs(y+x/2+47),1.5))-(x*(y-x+47)*cos(sqrt(abs(y-x+47))))/(2*pow(abs(y-x+47),1.5)));
  double dx = -((1-y/k)*(y/k+1)*(x/k+1)*((-y-47)*sin(sqrt(abs(x/2+y+47)))-x*sin(sqrt(abs(x-y-47)))))/k+((1-y/k)*(y/k+1)*(1-x/k)*((-y-47)*sin(sqrt(abs(x/2+y+47)))-x*sin(sqrt(abs(x-y-47)))))/k+(1-y/k)*(y/k+1)*(1-x/k)*(x/k+1)*(-sin(sqrt(abs(x-y-47)))-(x*(x-y-47)*cos(sqrt(abs(x-y-47))))/(2*pow(abs(x-y-47),1.5))+((-y-47)*cos(sqrt(abs(x/2+y+47)))*(x/2+y+47))/(4*pow(abs(x/2+y+47),1.5)));
  std::vector<double> grad;
  grad.push_back(dx);
  grad.push_back(dy);
  return grad;
}

double sum_abs(std::vector<double> gradient){
  return sqrt(pow(gradient[0],2)) + sqrt(pow(gradient[1],2));
}

// [[Rcpp::export]]
std::vector<double> gradientDescent(double starting_x, double starting_y, double gamma = 0.225, double epsilon = 0.001, int n = 5000, bool verbose = false) {
  double z_old = e(starting_x,starting_y);
  std::vector<double> par;
  par.push_back(starting_x);
  par.push_back(starting_y);
  for(int i =0; i <5000; i++){
    std::vector<double> gradient = e_prime(par);
    par[0] = par[0] + gamma * gradient[0];
    par[1] = par[1] + gamma * gradient[1];
    double z = e(par[0],par[1]);
    if(verbose){
      std::cout << "i = " << i << ", par = c(" << par[0] << par[1] << "), gradient = c(" << gradient[0] << gradient [1]<< "), size = " << sum_abs(gradient) << std::endl;
    }
    if(sum_abs(gradient) < epsilon || z_old > z){
      break;
    }
    z_old = z;
  }
  return par;
}