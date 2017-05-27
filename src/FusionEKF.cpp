#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;
  
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  
  Hj_ << 0, 0, 0, 0,
         0, 0, 0, 0,
        0, 0, 0, 0;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    
    previous_timestamp_ = measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = measurement_pack.raw_measurements_[0];
      double theta = measurement_pack.raw_measurements_[1];
      
      double p_x = rho * cos(theta);
      double p_y = rho * sin(theta);
      double v_x = 0;
      double v_y = 0;
      
      ekf_.x_ << p_x, p_y, v_x, v_y;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
             */
      double p_x = measurement_pack.raw_measurements_[0];
      double p_y = measurement_pack.raw_measurements_[1];
      ekf_.x_ << p_x, p_y, 0.0, 0.0;
    }
    
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1000, 0,
               0, 0, 0, 1000;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  double noise_ax = 9.0;
  double noise_ay = 9.0;
  
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
  
  
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;
  
  
  double dt_4_4 = (pow(dt, 4)/4);
  double dt_3_2 = (pow(dt, 3)/2);
  double dt_2 = pow(dt, 2);
  
  double dt_4_4_ax = dt_4_4 * noise_ax;
  double dt_4_4_ay = dt_4_4 * noise_ay;
  double dt_3_2_ax = dt_3_2 * noise_ax;
  double dt_3_2_ay = dt_3_2 * noise_ay;
  double dt_2_ax = dt_2 * noise_ax;
  double dt_2_ay = dt_2 * noise_ay;
  
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4_4_ax, 0, dt_3_2_ax, 0,
             0, dt_4_4_ay, 0, dt_3_2_ay,
             dt_3_2_ax, 0, dt_2_ax, 0,
             0, dt_3_2_ay, 0, dt_2_ay;
  
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.R_ = R_radar_;
    ekf_.H_ = Hj_;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    // Radar updates
  } else {
    // Laser updates
    ekf_.R_ = R_laser_;
    ekf_.H_ = MatrixXd(2, 4);
    ekf_.H_ << H_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
