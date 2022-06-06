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

/*  */
coupling_t forward_kinematics(UKF::Vector<3> position,
                                  UKF::Vector<3> orientation,
                                  calibration_t& calibration);

#endif