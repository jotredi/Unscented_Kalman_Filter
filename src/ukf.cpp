#include <iostream>
#include <cmath>
#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Wrap angle to pi
inline double wrapPI (double angle){
  // Args: angle in radians
  // Output: angle between (-pi, pi)
  while(angle < - M_PI) angle += 2*M_PI;
  while(angle > M_PI)   angle -= 2*M_PI;
  return angle;
}

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  //============================== NOISE ==============================//
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI/4;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //===================================================================//

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Augmented vector
  x_aug_ = VectorXd(n_aug_);

  // Augmented state covariance matrix
  P_aug_ = MatrixXd(n_aug_, n_aug_);

  // Sigma points matrix
  Xsig_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // Predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_/(lambda_ + n_aug_);
  for(int i=1; i<2*n_aug_+1; ++i) weights_(i) = 0.5/(n_aug_ + lambda_);

  //===================================================================//
  // Measurement matrix H - lidar
  H_lidar_ = MatrixXd(2, n_x_);
  H_lidar_.fill(0.0);
  H_lidar_(0,0) = H_lidar_(1,1) = 1;

  // Measurement noise covariance matrix R - lidar
  R_lidar_ = MatrixXd(2,2);
  R_lidar_.fill(0.0);
  R_lidar_(0,0) = pow(std_laspx_,2);
  R_lidar_(1,1) = pow(std_laspy_,2);

  // Measurement noise covariance matrix R - radar
  R_radar_ = MatrixXd(3,3);
  R_radar_.fill(0.0);
  R_radar_(0,0) = pow(std_radr_,2);
  R_radar_(1,1) = pow(std_radphi_,2);
  R_radar_(2,2) = pow(std_radrd_,2);

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  // Initialization
  if (!is_initialized_) {

    // Initialization of state x
    x_.fill(0.0);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

      // Recover ro & phi from measurement vector
      double ro = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];

      // Convert from polar to cartesian coordinates & initialize px, py
      x_(0) = ro * cos(phi);
      x_(1) = ro * sin(phi);

      // State Covariance Matrix
      P_.fill(0.0);
      // px & py initial variance
      P_(0,0) = P_(1,1) = 0.5;
      // v, yaw & yaw rate variance (high initial uncertainty)
      P_(2,2) = P_(3,3) = P_(4,4) = 1;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

      // Initialize state px & py from laser points
      x_(0) = meas_package.raw_measurements_[0];
      x_(1) = meas_package.raw_measurements_[1];

      // State Covariance Matrix
      P_.fill(0.0);
      // px & py initial variance
      P_(0,0) = P_(1,1) = 0.0225;
      // v, yaw & yaw rate variance (high initial uncertainty)
      P_(2,2) = P_(3,3) = P_(4,4) = 1;
    }

    // Set last timestamp
    previous_timestamp_ = meas_package.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************  PREDICTION  *****************************/

  // Calculate elapsed time & update last timestamp
  double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = meas_package.timestamp_;

  Prediction(dt);

  /*******************************  UPDATE  *******************************/

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    // Radar update
    UpdateRadar(meas_package);

  }else if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    // Laser update
    UpdateLidar(meas_package);
  }
}

void UKF::Prediction(double delta_t) {

  //================= GENERATE SIGMA POINTS =================//
  // Mean state
  x_aug_.fill(0.0);
  x_aug_.head(n_x_) = x_;

  // Covariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_,n_x_) = P_;
  P_aug_(5,5) = pow(std_a_,2);
  P_aug_(6,6) = pow(std_yawdd_,2);

  // Square root matrix
  MatrixXd L = P_aug_.llt().matrixL();

  // Create sigma points
  Xsig_.col(0) = x_aug_;
  for (int i=0; i<n_aug_; ++i){
    Xsig_.col(i+1) = x_aug_ + sqrt(n_aug_ + lambda_)*L.col(i);
    Xsig_.col(n_aug_+i+1) = x_aug_ - sqrt(n_aug_ + lambda_)*L.col(i);
  }

  //================= PREDICT SIGMA POINTS =================//
  // Noise vector
  VectorXd noise(n_x_);

  for (int j=0; j<Xsig_pred_.cols(); ++j){
    // Obtain values from Xsig
    double px = Xsig_(0,j);
    double py = Xsig_(1,j);
    double v = Xsig_(2,j);
    double yaw = Xsig_(3,j);
    double yaw_rate = Xsig_(4,j);
    double a_noise = Xsig_(5,j);
    double yawdd_noise = Xsig_(6,j);

    // Set noise vector
    noise(0) = pow(delta_t,2)*cos(yaw)*a_noise/2;
    noise(1) = pow(delta_t,2)*sin(yaw)*a_noise/2;
    noise(2) = delta_t*a_noise;
    noise(3) = pow(delta_t,2)*yawdd_noise/2;
    noise(4) = delta_t*yawdd_noise;

    // Check if yaw rate is zero
    if(yaw_rate){
      Xsig_pred_(0,j) = px + (v/yaw_rate)*(sin(yaw+yaw_rate*delta_t) - sin(yaw)) + noise(0);
      Xsig_pred_(1,j) = py + (v/yaw_rate)*(-cos(yaw+yaw_rate*delta_t) + cos(yaw)) + noise(1);
      Xsig_pred_(2,j) = v + noise(2);
      Xsig_pred_(3,j) = yaw + yaw_rate*delta_t + noise(3);
      Xsig_pred_(4,j) = yaw_rate + noise(4);
    }
    else{
      Xsig_pred_(0,j) = px + v*cos(yaw)*delta_t + noise(0);
      Xsig_pred_(1,j) = py + v*sin(yaw)*delta_t + noise(1);
      Xsig_pred_(2,j) = v + noise(2);
      Xsig_pred_(3,j) = yaw + noise(3);
      Xsig_pred_(4,j) = yaw_rate + noise(4);
    }
  }

  //================= PREDICT MEAN & COVARIANCE =================//
  // Predict state mean
  x_.fill(0.0);
  for (int j=0; j<Xsig_pred_.cols(); ++j){
    x_ += weights_(j) * Xsig_pred_.col(j);
  }

  // Predict state covariance matrix
  P_.fill(0.0);
  for (int j=0; j<Xsig_pred_.cols(); ++j){
    // state difference
    VectorXd diff = Xsig_pred_.col(j) - x_;
    // angle normalize
    diff(3) = wrapPI(diff(3));

    P_ += weights_(j) * diff * diff.transpose();
  }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {

  // Measurement vector
  VectorXd z = meas_package.raw_measurements_;

  // Calculate error vector y
  VectorXd y = z - H_lidar_ * x_;

  // Matrix S
  MatrixXd S = H_lidar_ * P_ * H_lidar_.transpose() + R_lidar_;

  // Kalman gain
  MatrixXd K = P_ * H_lidar_.transpose() * S.inverse();

  // State update
  x_ += K * y;

  // Covariance matrix update
  MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
  P_ = (I - K * H_lidar_) * P_;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {

  //================= PREDICT RADAR MEASUREMENT =================//
  // Set measurement dimension
  int n_z = 3;

  // Matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  // Measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  // Transform sigma points into measurement space
  for (int j=0; j<Xsig_pred_.cols(); ++j){
    // Get state values
    double px = Xsig_pred_(0,j);
    double py = Xsig_pred_(1,j);
    double v = Xsig_pred_(2,j);
    double yaw = Xsig_pred_(3,j);

    // Set new sigma points
    Zsig(0,j) = sqrt(pow(px,2) + pow(py,2));
    Zsig(1,j) = atan2(py,px);
    Zsig(2,j) = (px*cos(yaw)*v + py*sin(yaw)*v)/Zsig(0,j);
  }

  // Calculate mean predicted measurement
  z_pred.fill(0.0);
  for(int j=0; j<Zsig.cols(); ++j){
    z_pred += weights_(j) * Zsig.col(j);
  }

  // Calculate innovation covariance matrix S
  S.fill(0.0);
  for(int j=0; j<Zsig.cols(); ++j){
    // meas difference
    VectorXd diff = Zsig.col(j) - z_pred;

    // angle normalize
    diff(1) = wrapPI(diff(1));

    S += weights_(j) * diff * diff.transpose();
  }

  // Add additive noise
  S += R_radar_;

  //====================== UPDATE ======================//
  // Measurement vector
  VectorXd z = meas_package.raw_measurements_;

  // Cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int j=0; j<Xsig_pred_.cols(); ++j){
    // State difference
    VectorXd x_diff = Xsig_pred_.col(j) - x_;
    // angle normalize
    x_diff(3) = wrapPI(x_diff(3));

    // Meas difference
    VectorXd z_diff = Zsig.col(j) - z_pred;
    // angle normalize
    z_diff(1) = wrapPI(z_diff(1));

    Tc += weights_(j) * x_diff * z_diff.transpose();
  }

  // Kalman gain
  MatrixXd K = Tc * S.inverse();

  // Measurement error vector
  VectorXd y = z - z_pred;
  // angle normalize
  y(1) = wrapPI(y(1));

  // Update state mean & covariance matrix
  x_ += K * y;
  P_ -= K * S * K.transpose();
}
