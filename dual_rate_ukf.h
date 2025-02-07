#ifndef INTERFACE_H
#define INTERFACE_H

/* Calibration struct definition */
typedef struct pre_calibration {
    UKF::Matrix<3, 3> d_source_pos;
    UKF::Matrix<3, 3> d_source_moment;
    UKF::Matrix<3, 3> d_sensor_pos;
    UKF::Matrix<3, 3> d_sensor_moment;
} pre_calibration_t;

typedef struct calibration {
    UKF::Matrix<3, 1> d_source_pos[3];
    UKF::Matrix<3, 1> d_source_moment[3];
    UKF::Matrix<3, 1> d_sensor_pos[3];
    UKF::Matrix<3, 1> d_sensor_moment[3];
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

Eigen::Matrix<real_t, Eigen::Dynamic, 6> run_high_ukf(
    Eigen::Matrix<real_t, Eigen::Dynamic, 9> &coupling_hi);

#endif