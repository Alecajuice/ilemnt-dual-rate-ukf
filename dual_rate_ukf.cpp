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
#include <string>

calibration_t calibration_hi;
calibration_t calibration_lo;

Eigen::Affine3d pose2trans(UKF::Vector<3> position, UKF::Vector<3> orientation) {
    real_t mag = orientation.norm();
    if (mag == 0) {
        // TODO: what to do in this case
        return Eigen::Translation3d(position) * Eigen::AngleAxisd(0, Eigen::Vector3d(1, 0, 0));
    }
    return Eigen::Translation3d(position) * Eigen::AngleAxisd(mag, orientation / mag);
}

coupling_t forward_kinematics(UKF::Vector<3> position,
                                  UKF::Vector<3> orientation,
                                  calibration_t& calibration)
{
    coupling_t coupling = coupling_t::Zero();
    Eigen::Affine3d P = pose2trans(position, orientation);
    Eigen::Affine3d R = Eigen::Translation3d(-position) * P; // zero out translation portion of P
    for (int j = 0; j < 3; j++) { // sensor index
        // sensor position and moment in source coordinates
        Eigen::Vector3d sensor_pos = P * calibration.d_sensor_pos(Eigen::all, j);
        Eigen::Vector3d sensor_moment = R * calibration.d_sensor_moment(Eigen::all, j);
        for (int i = 0; i < 3; i++) { // source index

            // vector distance r between source and sensor in source coordinates
            Eigen::Vector3d rSoSe = sensor_pos - calibration.d_source_pos(Eigen::all, i);
            // vector r magnitude
            real_t rSoSeMag = rSoSe.norm() == 0 ? 1e-10 : rSoSe.norm();
            // TODO: what to do in this case
            // vector r unit
            Eigen::Vector3d rSoSeUnit = rSoSe / rSoSeMag;

            // Find magnetic field vector at the sensor location, using the dipole model
            Eigen::Vector3d B = (1 / std::pow(rSoSeMag, 3)) *
                (3 * calibration.d_source_moment(Eigen::all, i).dot(rSoSeUnit) * rSoSeUnit -
                    calibration.d_source_moment(Eigen::all, i));

            // coupling matrix using antenna model
            coupling(3 * i + j) = B.dot(sensor_moment) / 250;
        }
    }
    // std::cout << coupling << std::endl;
    return coupling;
}

int main() {
    std::cout << "Starting dual UKF..." << std::endl;
    // Inputs
    Eigen::Matrix<real_t, Eigen::Dynamic, 9> couplings_hi_at_lo;
    Eigen::Matrix<real_t, Eigen::Dynamic, 9> couplings_lo;
    Eigen::Matrix<real_t, Eigen::Dynamic, 9> couplings_hi;
    // Outputs
    Eigen::Matrix<real_t, Eigen::Dynamic, 9> biases;
    Eigen::Matrix<real_t, Eigen::Dynamic, 6> poses;
    
    std::string trace = "aluminum_sheet";
    try {
        // Read calibration files
        matio::read_mat("XYZ_hr_cal.mat", "d_source_pos", calibration_hi.d_source_pos);
        matio::read_mat("XYZ_hr_cal.mat", "d_source_moment", calibration_hi.d_source_moment);
        matio::read_mat("XYZ_hr_cal.mat", "d_sensor_pos", calibration_hi.d_sensor_pos);
        matio::read_mat("XYZ_hr_cal.mat", "d_sensor_moment", calibration_hi.d_sensor_moment);
        matio::read_mat("XYZ_lr_cal.mat", "d_source_pos", calibration_lo.d_source_pos);
        matio::read_mat("XYZ_lr_cal.mat", "d_source_moment", calibration_lo.d_source_moment);
        matio::read_mat("XYZ_lr_cal.mat", "d_sensor_pos", calibration_lo.d_sensor_pos);
        matio::read_mat("XYZ_lr_cal.mat", "d_sensor_moment", calibration_lo.d_sensor_moment);

        // Read coupling files
        matio::read_mat(trace + "_couplings.mat", "couplings_hi_at_lo", couplings_hi_at_lo);
        matio::read_mat(trace + "_couplings.mat", "couplings_lo", couplings_lo);
        matio::read_mat(trace + "_couplings.mat", "couplings_hi", couplings_hi);

        // Read bias file
        matio::read_mat(trace + "_biases.mat", "biases", biases);
    }
    catch (const std::exception & ex) {
        std::cout << "error:" << ex.what() << std::endl;
    }

    // UKF
    std::cout << "Running low rate UKF..." << std::endl;
    // biases = run_low_ukf(couplings_hi_at_lo, couplings_lo);

    Eigen::Matrix<real_t, Eigen::Dynamic, 9> unbiased_hi =
        couplings_hi -
            biases(Eigen::VectorXi::LinSpaced(128 * biases.rows(), 0, biases.rows() - 1), Eigen::all);

    std::cout << "Running high rate UKF..." << std::endl;
    poses = run_high_ukf(unbiased_hi);
    try {
        std::cout << "Writing biases to file..." << std::endl;
        matio::write_mat(trace + "_biases.mat", "biases", biases, true);
        std::cout << "Writing poses to file..." << std::endl;
        matio::write_mat(trace + "_posesukf_no_dynamics.mat", "posesukf", poses, true);
    }
    catch (const std::exception & ex) {
        std::cout << "error:" << ex.what() << std::endl;
    }

    return 0;
}