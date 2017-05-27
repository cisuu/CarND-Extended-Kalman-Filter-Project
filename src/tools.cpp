#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  
  int size = estimations.size();
 
  if(size != ground_truth.size() || size == 0){
    return rmse;
  }
  
  for(int i=0; i < size; ++i){
    // ... your code here
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }
  rmse = rmse/size;
  
  rmse = rmse.array().sqrt();
  
  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);
  
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);
  
  if(fabs(px) < 0.001 || fabs(py) < 0.001){
    cout << "error" << endl;
    return Hj;
  }
  
  double px_py_pow_sum = pow(px,2) + pow(py,2);
  double Hj_0_0 = px/(sqrt(px_py_pow_sum));
  double Hj_0_1 = py/(sqrt(px_py_pow_sum));
  double Hj_0_2 = 0;
  double Hj_0_3 = 0;
  double Hj_1_0 = -(py/px_py_pow_sum);
  double Hj_1_1 = px/px_py_pow_sum;
  double Hj_1_2 = 0;
  double Hj_1_3 = 0;
  double Hj_2_0 = (py*(vx*py-vy*px))/pow(px_py_pow_sum, 1.5);
  double Hj_2_1 = (px*(vy*px-vx*py))/pow(px_py_pow_sum, 1.5);
  double Hj_2_2 = px/(sqrt(px_py_pow_sum));
  double Hj_2_3 = py/(sqrt(px_py_pow_sum));
  
  if(fabs(px_py_pow_sum) < 0.001){
    cout << "CalculateJacobian () - Error - Division by Zero" << endl;
    return Hj;
  }
  
  //compute the Jacobian matrix
  
  Hj << Hj_0_0, Hj_0_1, Hj_0_2, Hj_0_3,
        Hj_1_0, Hj_1_1, Hj_1_2, Hj_1_3,
        Hj_2_0, Hj_2_1, Hj_2_2, Hj_2_3;
  
  cout << "Hj = " << Hj << endl;
  
  return Hj;
}
