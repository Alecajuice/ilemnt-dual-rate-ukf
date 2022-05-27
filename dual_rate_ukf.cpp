#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <cmath>
#define UKF_DOUBLE_PRECISION
#include "UKF/Types.h"
#include "UKF/Integrator.h"
#include "UKF/StateVector.h"
#include "UKF/MeasurementVector.h"
#include "UKF/Core.h"
#include "dual_rate_ukf.h"
#include "MATio"
#include <iostream>

calibration_t calibration_hi;
calibration_t calibration_lo;

coupling_t forward_kinematics(UKF::Vector<3> position,
                                  UKF::Vector<3> orientation,
                                  calibration_t calibration)
{
    return coupling_t::Zero();
}

int main() {
    std::cout << "Starting dual UKF..." << std::endl;
    Eigen::Matrix3f ff;
  
    // read matrix 'ff' from file 'data.mat' in as floats
    try {
        matio::read_mat("data.mat", "test", ff);
        std::cout << "test=" << ff << std::endl;
    
        // write it back out as double precision 'dd'
        matio::write_mat("data.mat", "dd", ff.cast<double>());
    }
    catch (const std::exception & ex) {
        std::cout << "error:" << ex.what() << std::endl;
    }
    return 0;
}