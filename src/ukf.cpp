#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // state dim
  n_x_ = 5;
  // augment. state dim, add 2 noise parameters
  n_aug_ = n_x_ + 2;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // For a Gaussian distribution, we expect the acceleration to be between
  // −60​ m/s2​​ and +60​m/s2 or −60​rad/s2​​ and +60​rad/s2​​ ninety-five percent of the time .
  // Process noise standard deviation longitudinal acceleration in m/s^2
  // It means 60 m/s2 = 2*std_a -> std_a = 60/2 = 30
  // the fastest measured linear acceleration for a street legal sports car
  // is currently 0 to 60 mph in 2.2 seconds. 0 to 60 mph in 2.2 seconds is about 12​m/s2
  // 12m/s2 = 2*std_a -> std_a = 12/2 = 6
  std_a_ = 10;    //given was 30

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 4;  //given was 30

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;

  time_us_ = 0;

  Xsig_pred_= MatrixXd(n_x_, 2 * n_aug_ + 1);

  lambda_= 3 - n_aug_;

  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i=1; i<2*n_aug_ + 1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  H_laser_ = MatrixXd(2, 5);
  R_laser_ = MatrixXd(2,2);
  H_laser_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0;
  R_laser_ << 0.0225, 0,
          0, 0.0225;

  NIS_radar_ = 0;

  NIS_laser_ = 0;
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
  if (!is_initialized_) {
    // initialize state vector x and covariance matrix P
    P_ = MatrixXd::Identity(n_x_, n_x_);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      cout << "\n Incoming RADAR\n----------------" << endl;
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];

      double px = rho * cos(phi);
      double py = rho * sin(phi);
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);
      double v = sqrt(vx*vx + vy*vy);
      double psi = atan2(vy, vx);   // or 0?

      x_ << px, py, v, psi, 0;

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      cout << "\n Incoming LASER\n----------------" << meas_package.sensor_type_<<"\n"<< endl;
      cout << meas_package.raw_measurements_ << endl;


      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }
    else
    {
      cout<<"WARNING: unknown sensor type!"<<endl;
    }
    time_us_ = meas_package.timestamp_ ;
    is_initialized_ = true;
    cout << "x_: \n"<< x_ <<"\n"<< endl;
    return;
  }

  double dt = (meas_package.timestamp_ - time_us_) / 1000000.;
  time_us_ = meas_package.timestamp_;
  if (dt > 0.05)
    cout << "### Caution, dt is: "<< dt <<"\n"<< endl;

  Prediction(dt);

  if(meas_package.sensor_type_ == MeasurementPackage::RADAR ){

    if(use_radar_){

      UpdateRadar(meas_package);
    }
    else
      cout<<"WARNING: Ignoring radar!"<<endl;
  } else if(meas_package.sensor_type_ == MeasurementPackage::LASER )
  {

    if(use_laser_){

      UpdateLidar(meas_package);
    }
    else
      cout<<"WARNING: Ignoring laser!"<<endl;
  }else
  {
    cout<<"WARNING: unknown sensor type!"<<endl;
  }

  cout << "x_ "<< x_ <<"\n"<< endl;
  cout << "P_ "<< P_ <<"\n"<< endl;
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
  // 1. Generate augmented Sigma points
  //--------------------------------------------
  //create augmented mean vector
  static VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  static MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  static MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  static MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+ n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+ n_aug_) * L.col(i);
  }

  //----------------------------------
  // 2. Predict augmented Sigma points
  //----------------------------------
  double p_x, p_y, v, yaw, yawd, nu_a, nu_yawdd, px_p, py_p;
  //predict sigma points
  for (int i = 0; i< 2* n_aug_ +1; i++)
  {
    //extract values for better readability
    //YW: from Xsig_aug, not from x_aug
    p_x = Xsig_aug(0,i);
    p_y = Xsig_aug(1,i);
    v = Xsig_aug(2,i);
    yaw = Xsig_aug(3,i);
    yawd = Xsig_aug(4,i);
    nu_a = Xsig_aug(5,i);
    nu_yawdd = Xsig_aug(6,i);


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

  //------------------------------------
  // 3. Predict mean and covariance matrix
  //------------------------------------
  static VectorXd x_pred = VectorXd(n_x_);
  static MatrixXd P_pred = MatrixXd(n_x_, n_x_);
  x_pred.fill(0.);
  P_pred.fill(0.);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_pred = x_pred + weights_(i) * Xsig_pred_.col(i);
  }

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_pred = P_pred + weights_(i) * x_diff * x_diff.transpose() ;
  }
  x_ = x_pred;
  P_ = P_pred;
}

/*void UKF::PredictionLidar(double delta_t){
  double noise_ax, noise_ay = 9.;

  cout << "\n LIDAR prediction \n----------------" << endl;
  F_ << 1, 0, delta_t, 0,
         0, 1, 0,  delta_t,
         0, 0, 1,  0,
         0, 0, 0,  1;
  // Set the covariance matrix
  Q_ << pow(delta_t,4)/4*noise_ax, 0, pow(delta_t,3)/2*noise_ax, 0,
         0, pow(delta_t,4)/4*noise_ay, 0, pow(delta_t,3)/2*noise_ay,
         pow(delta_t,3)/2*noise_ax, 0, pow(delta_t,2)*noise_ax, 0,
         0, pow(delta_t,3)/2*noise_ay, 0, pow(delta_t,2)*noise_ay;
  cout << "\n LIDAR prediction 2\n----------------" << endl;
  VectorXd x = VectorXd(4);
  x(0) << x_(0), x(1), 0, 0;

  x = F_ * x;
  cout << "\n LIDAR prediction 3\n----------------" << endl;
  P_ = F_ * P_ * F_.transpose() + Q_;
}*/
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
  cout << "\n LIDAR update \n----------------" << endl;

  VectorXd z_pred = H_laser_ * x_;

  VectorXd y = meas_package.raw_measurements_ - z_pred;

  MatrixXd Ht = H_laser_.transpose();
  MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;

  //new state
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);

  P_ = (I - K * H_laser_) * P_;
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
  // -------------------
  // Predict measurement
  // -------------------
  cout << "\n RADAR update \n----------------" << endl;
  //create matrix for sigma points in measurement space, 3 rows for radar dim
  int n_z = 3;
  static MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  =  Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;
    // avoid zero devision
    if(fabs(p_x) <= 0.0001){
            p_x = 0.0;
    }
    if(fabs(p_y) <= 0.0001){
            p_y = 0.0;
    }
    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  static VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_ + 1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  static MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  static MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_ * std_radr_, 0, 0,
          0, std_radphi_ * std_radphi_, 0,
          0, 0,std_radrd_ * std_radrd_;
  S = S + R;

  // -----------------
  // Update State
  // -----------------
  //create matrix for cross correlation Tc
  static MatrixXd Tc = MatrixXd(n_x_, n_z);
  //calculate cross correlation matrix
  Tc.fill(0.);
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

  //Kalman gain K;
  static MatrixXd K = MatrixXd(n_x_, n_z);
  K = Tc * S.inverse();

  //residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}
