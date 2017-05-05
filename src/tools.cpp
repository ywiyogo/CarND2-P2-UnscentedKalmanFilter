#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
    VectorXd rmse = VectorXd(4);
    rmse << 0, 0, 0, 0;
    // Validity check
    if(estimations.size() == 0 || estimations.size() != ground_truth.size()){
        cout << "Invalid data" << endl;
        return rmse;
    }

    //accumulate squared residuals
    for(unsigned int i=0; i < estimations.size(); ++i){
        VectorXd diff = estimations[i] - ground_truth[i];
        diff = diff.array()*diff.array();
        rmse += diff;
    }

    //calculate mean
    //cout << "estimations.size(): "<< estimations.size()<<endl;
    rmse /= estimations.size();
    //squared root
    rmse = rmse.array().sqrt();
    return rmse;
}
