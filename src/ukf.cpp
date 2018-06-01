#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 6;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.8;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:
  
  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_=false;
  n_x_=5;
  n_aug_=7;
  lambda_=3-n_aug_;
  Xsig_pred_=MatrixXd(n_x_,2*n_aug_+1);
  weights_=VectorXd(2*n_aug_+1);
  
  P_<<0.15,0,0,0,0,
      0,0.15,0,0,0,
      0,0,10,0,0,
      0,0,0,10,0,
      0,0,0,0,20;
  
  
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if(!is_initialized_){
     if(meas_package.sensor_type_==MeasurementPackage::LASER)
     {
       x_<<meas_package.raw_measurements_[0],meas_package.raw_measurements_[1],1,0,.5;
       time_us_=meas_package.timestamp_;
      }
    else
    {
      double pos_x=meas_package.raw_measurements_[0]*cos(meas_package.raw_measurements_[1]*180/3.14);
      double pos_y=meas_package.raw_measurements_[0]*sin(meas_package.raw_measurements_[1]*180/3.14);
      x_<<pos_x,pos_y,0,0,0;
      time_us_=meas_package.timestamp_;
    }
     is_initialized_ = true;
    // cout<<"X=\n"<<x_;
   // cout<<"P=\n"<<P_;
      return;
  }
  
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
	time_us_ = meas_package.timestamp_;
  //cout<<"DT="<<dt<<"\n";
  //cout<<"Predicting\n";
  Prediction(dt);
  //cout<<"END Prediction\n";
  if(meas_package.sensor_type_==MeasurementPackage::RADAR)
  {
    if(use_radar_){
      //cout<<"Radar Update\n";
      UpdateRadar(meas_package);
     // cout<<"Radar Updation Completed\n";   
    }
  }
  else if(meas_package.sensor_type_==MeasurementPackage::LASER){
    if(use_laser_){
      //cout<<"LIDAR Update\n";
      UpdateLidar(meas_package);
      //cout<<"LIDAR Updation Completed\n";
    }
  }
  //cout<<"X=\n"<<x_;
//cout<<"P=\n"<<P_;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
    
    //Augmented State
    
    VectorXd x_aug=VectorXd(n_aug_);
    x_aug.head(5)=x_;
    x_aug(5)=0.0;
    x_aug(6)=0.0;
    
    Xsig_pred_.fill(0.0);
    //Augmented Covariance Matrix
    MatrixXd P_aug=MatrixXd(n_aug_,n_aug_);
    P_aug.fill(0.0);
    P_aug.topLeftCorner(5,5)=P_;
    P_aug(5,5)=std_a_*std_a_;
    P_aug(6,6)=std_yawdd_*std_yawdd_;
     

    MatrixXd L = P_aug.llt().matrixL();
    MatrixXd Xsig_aug=MatrixXd(n_aug_,2*n_aug_+1);
    Xsig_aug.fill(0.0);
    Xsig_aug.col(0)=x_aug;
    for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
  //cout<<"XSig_aug=\n"<<Xsig_aug<<"\n";
  
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
 
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights_
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }
  
  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_+ weights_(i) * Xsig_pred_.col(i);
  }
   
  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
  
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  
  VectorXd z=VectorXd(2);
   z<<meas_package.raw_measurements_[0],meas_package.raw_measurements_[1];
  
  MatrixXd H_=MatrixXd(2,5);
  MatrixXd R_=MatrixXd(2,2);
  R_<<std_laspy_*std_laspy_,0,
      0,std_laspx_*std_laspx_;


  H_<<1,0,0,0,0,
      0,1,0,0,0;
  VectorXd z_pred = H_ * x_;
  
	VectorXd y = z - z_pred;
  
	MatrixXd Ht = H_.transpose();
  
	MatrixXd S = H_ * P_ * Ht + R_;
  
	MatrixXd Si = S.inverse();
  
	MatrixXd PHt = P_ * Ht;

	MatrixXd K = PHt * Si;
  
	//new estimate

	x_ = x_ + (K * y);

	long x_size = x_.size();

	MatrixXd I = MatrixXd::Identity(x_size, x_size);

	P_ = (I - K * H_) * P_;
  VectorXd NIS=y.transpose()*Si*y;
  cout<<"LIDAR NIS="<<NIS<<",\n";
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  //cout<<"Xsig_pred_="<<Xsig_pred_<<"\n";
  int n_z=3;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);


  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }
  //cout<<"Zsig="<<Zsig<<"\n";
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2 * n_aug_ + 1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  //cout<<"weights="<<weights_<<"\n";
  //cout<<"z_pred="<<z_pred<<"\n";
  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  //cout<<"S1="<<S<<"\n";
  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0.0, 0.0,
          0.0, std_radphi_*std_radphi_, 0.0,
          0.0, 0.0,std_radrd_*std_radrd_;
  S = S + R;
  
  //cout<<"S="<<S<<"\n";
  
 
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  //cout<<"Tc="<<Tc<<"\n";
  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  //cout<<"K="<<K<<"\n";
  //residual
  VectorXd z=VectorXd(3);
  z<<meas_package.raw_measurements_[0],meas_package.raw_measurements_[1],meas_package.raw_measurements_[2];
  
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  
  x_ = x_ + K * z_diff;
  //cout<<"X="<<x_<<"\n";
  
  P_ = P_ - K*S*K.transpose();
  //cout<<"P="<<P_<<"\n";
  VectorXd NIS=z_diff.transpose()*S.inverse()*z_diff;
  cout<<"RADAR NIS="<<NIS<<",\n";
}
