#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) 
{
  /**
    * Calculate the RMSE here.
  */
   	VectorXd rmse(4);
	VectorXd xest(4), xtrue(4), xdiff(4), xdiff2(4);
	rmse << 0,0,0,0;

    int n = estimations.size();
	if(estimations.size() == 0 || estimations.size() != ground_truth.size())
    {
        std::cout << "Invalid estimation or ground truth" << std::endl;
	    return rmse;
    }

	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
        // ... your code here
		xest = estimations[i];
		xtrue = ground_truth[i];
		xdiff = (xest - xtrue);
		xdiff2 = xdiff.array()*xdiff.array();
		//rmse = rmse +  xdiff.array() * xdiff.array();

		rmse += xdiff2;
	}

	//calculate the mean
	rmse = rmse *1.0/(float)n;

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

  float px2py2 = px*px + py*py;

	//check division by zero
	if(px2py2 == 0.0)
	{
	  cout << "Error" << endl;
	}
	else
	//compute the Jacobian matrix
  {
    float sqrtpx2py2 = sqrt(px2py2);
    float px2py2_32 = pow(px2py2, 1.5);
    Hj << px/sqrtpx2py2 , py/sqrtpx2py2, 0.0, 0.0,
          -py/px2py2, px/px2py2, 0.0, 0.0,
          py*(vx*py - vy*px)/px2py2_32,
          px*(vy*px - vx*py)/px2py2_32,
          px/sqrtpx2py2, py/sqrtpx2py2;
  }
	return Hj;
}
