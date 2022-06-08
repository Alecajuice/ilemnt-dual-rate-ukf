#ifndef INTERFACE_H
#define INTERFACE_H

/* Calibration struct definition */
typedef struct calibration {
    UKF::Matrix<3, 3> d_source_pos;
    UKF::Matrix<3, 3> d_source_moment;
    UKF::Matrix<3, 3> d_sensor_pos;
    UKF::Matrix<3, 3> d_sensor_moment;
} calibration_t;

extern calibration_t calibration_hi;
extern calibration_t calibration_lo;

/* Coupling type */
typedef UKF::Vector<9> coupling_t;

/* FK */
coupling_t forward_kinematics(UKF::Vector<3> position,
                                  UKF::Vector<3> orientation,
                                  calibration_t& calibration);

/* Run UKF functions */
Eigen::Matrix<real_t, Eigen::Dynamic, 9> run_low_ukf(
    Eigen::Matrix<real_t, Eigen::Dynamic, 9> &coupling_hi_at_lo,
    Eigen::Matrix<real_t, Eigen::Dynamic, 9> &coupling_lo);

#endif