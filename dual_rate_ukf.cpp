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

coupling_t forward_kinematics(UKF::Vector<3> position,
                                  UKF::Vector<3> orientation,
                                  calibration_t& calibration)
{
    coupling_t coupling = coupling_t::Zero();
    real_t mag = orientation.norm();
    Eigen::Affine3d R(Eigen::AngleAxisd(mag, orientation / mag)); // rotation matrix w/ only orientation
    Eigen::Affine3d P = Eigen::Translation3d(position) * R; // transform w/ both pos and orientation
    for (int j = 0; j < 3; j++) { // sensor index
        // sensor position and moment in source coordinates
        Eigen::Vector3d sensor_pos = P * calibration.d_sensor_pos[j];
        Eigen::Vector3d sensor_moment = R * calibration.d_sensor_moment[j];
        for (int i = 0; i < 3; i++) { // source index
            // vector distance r between source and sensor in source coordinates
            Eigen::Vector3d rSoSe = sensor_pos - calibration.d_source_pos[i];
            // vector r magnitude
            real_t rSoSeMag = rSoSe.norm();
            // vector r unit
            Eigen::Vector3d rSoSeUnit = rSoSe / rSoSeMag;

            // Find magnetic field vector at the sensor location, using the dipole model
            auto so_mo = calibration.d_source_moment[i];
            Eigen::Vector3d B = (1 / std::pow(rSoSeMag, 3)) *
                (3 * so_mo.dot(rSoSeUnit) * rSoSeUnit - so_mo);

            // coupling matrix using antenna model
            coupling(3 * i + j) = B.dot(sensor_moment) / 250;
        }
    }
    // std::cout << coupling << std::endl;
    return coupling;
}

void preprocess_calibration(pre_calibration_t &pre, calibration_t &cal) {
    for (int i = 0; i < 3; i++) {
        cal.d_source_pos[i] = pre.d_source_pos(Eigen::all, i);
        cal.d_source_moment[i] = pre.d_source_moment(Eigen::all, i);
        cal.d_sensor_pos[i] = pre.d_sensor_pos(Eigen::all, i);
        cal.d_sensor_moment[i] = pre.d_sensor_moment(Eigen::all, i);
    }
}

int main() {
    std::cout << "Starting dual UKF..." << std::endl;
    // Inputs
    Eigen::Matrix<real_t, Eigen::Dynamic, 9> couplings_hi_at_lo;
    Eigen::Matrix<real_t, Eigen::Dynamic, 9> couplings_lo;
    Eigen::Matrix<real_t, Eigen::Dynamic, 9> couplings_hi;
    pre_calibration_t pre_hi, pre_lo;
    // Outputs
    Eigen::Matrix<real_t, Eigen::Dynamic, 9> biases;
    Eigen::Matrix<real_t, Eigen::Dynamic, 6> poses;
    
    std::string trace = "aluminum_sheet";
    try {
        // Read calibration files
        matio::read_mat("XYZ_hr_cal.mat", "d_source_pos", pre_hi.d_source_pos);
        matio::read_mat("XYZ_hr_cal.mat", "d_source_moment", pre_hi.d_source_moment);
        matio::read_mat("XYZ_hr_cal.mat", "d_sensor_pos", pre_hi.d_sensor_pos);
        matio::read_mat("XYZ_hr_cal.mat", "d_sensor_moment", pre_hi.d_sensor_moment);
        matio::read_mat("XYZ_lr_cal.mat", "d_source_pos", pre_lo.d_source_pos);
        matio::read_mat("XYZ_lr_cal.mat", "d_source_moment", pre_lo.d_source_moment);
        matio::read_mat("XYZ_lr_cal.mat", "d_sensor_pos", pre_lo.d_sensor_pos);
        matio::read_mat("XYZ_lr_cal.mat", "d_sensor_moment", pre_lo.d_sensor_moment);

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

    // Preprocess calibrations
    preprocess_calibration(pre_hi, calibration_hi);
    preprocess_calibration(pre_lo, calibration_lo);

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