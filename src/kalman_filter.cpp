#include "kalman_filter.h"
#include <iostream>
#include <cmath>

#define EPS 0.001

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd y = z - H_ * x_;
  StateCovarianceUpdate(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  VectorXd h = VectorXd(3);
  double p_x = x_(0);
  double p_y = x_(1);
  double v_x = x_(2);
  double v_y = x_(3);
  double rho = sqrt(pow(p_x, 2) + pow(p_y, 2));
  if(fabs(p_x) < EPS) {
    p_x = EPS;
  }
  double theta = atan2(p_y, p_x);
  if (fabs(rho) < EPS){
    rho = EPS;
  }
  double rho_dot = (p_x * v_x + p_y * v_y) / (rho);

  h << rho, theta, rho_dot;
  
  VectorXd y = z - h;
  if(fabs(y[1]) > M_PI){
      y[1] = atan2(sin(y[1]), cos(y[1]));
  }
  StateCovarianceUpdate(y);
}

void KalmanFilter::StateCovarianceUpdate(const VectorXd &y){
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;
  
  x_ = x_ + (K* y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}


