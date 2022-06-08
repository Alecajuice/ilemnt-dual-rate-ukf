#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#define UKF_DOUBLE_PRECISION
#include "UKF/Types.h"
#include "UKF/Integrator.h"
#include "UKF/StateVector.h"
#include "UKF/MeasurementVector.h"
#include "UKF/Core.h"
#include "dual_rate_ukf.h"
#include <iostream>

/* State vector definition */
enum Low_States {
    Position,
    Orientation,
    Coupling_Bias
};

using Low_StateVector = UKF::StateVector<
    UKF::Field<Position, UKF::Vector<3>>,
    UKF::Field<Orientation, UKF::Vector<3>>,
    UKF::Field<Coupling_Bias, coupling_t>
>;

/* State transition function */
namespace UKF {

/* SFWA state vector process model. */
template <> template <>
Low_StateVector Low_StateVector::derivative<>() const {
    Low_StateVector output;
    output.set_field<Position>(UKF::Vector<3>(0, 0, 0));
    output.set_field<Orientation>(UKF::Vector<3>(0, 0, 0));
    output.set_field<Coupling_Bias>(coupling_t::Zero());
    return output;
}

}

/* Measurement vector definition */
enum Low_Measurements {
    High_At_Low_Coupling,
    Low_Coupling
};

using Low_MeasurementVector = UKF::DynamicMeasurementVector<
    UKF::Field<High_At_Low_Coupling, coupling_t>,
    UKF::Field<Low_Coupling, coupling_t>
>;

/* Measurement functions */
namespace UKF {

template <> template <>
coupling_t Low_MeasurementVector::expected_measurement
<Low_StateVector, High_At_Low_Coupling>(const Low_StateVector& state) {
    return forward_kinematics(state.get_field<Position>(),
                              state.get_field<Orientation>(),
                              calibration_hi)
            + state.get_field<Coupling_Bias>();
}

template <> template <>
coupling_t Low_MeasurementVector::expected_measurement
<Low_StateVector, Low_Coupling>(const Low_StateVector& state) {
    return forward_kinematics(state.get_field<Position>(),
                              state.get_field<Orientation>(),
                              calibration_lo);
}

}

/* UKF Core definition */
using Low_Core = UKF::Core<
    Low_StateVector,
    Low_MeasurementVector,
    UKF::IntegratorRK4
>;

/* Run UKF */
Eigen::Matrix<real_t, Eigen::Dynamic, 9> run_low_ukf(
    Eigen::Matrix<real_t, Eigen::Dynamic, 9> &coupling_hi_at_lo,
    Eigen::Matrix<real_t, Eigen::Dynamic, 9> &coupling_lo)
{
    // Return value
    Eigen::Matrix<real_t, Eigen::Dynamic, 9> biases(coupling_lo.rows(), 9);

    // Const values
    const real_t initial_variance = 1e-6;

    const real_t base_process_noise = 2e-3;
    const real_t process_noise_trans = base_process_noise;
    const real_t process_noise_rot = base_process_noise * 2;
    const real_t process_noise_bias = 1e-5;

    const real_t measurement_noise_lo = 1e-7;
    const real_t measurement_noise_hi = measurement_noise_lo * 100;

    const real_t dt = 1.0/1500;

    // Initialize filter
    Low_Core filter;

    filter.state.set_field<Position>(UKF::Vector<3>::Zero());
    filter.state.set_field<Orientation>(UKF::Vector<3>::Zero());
    filter.state.set_field<Coupling_Bias>(UKF::Vector<9>::Zero());
    std::cout << filter.state.get_field<Coupling_Bias>() << std::endl;

    filter.covariance = Low_StateVector::CovarianceMatrix::Identity() * initial_variance;
    std::cout << filter.covariance << std::endl;

    Eigen::Vector<real_t, 15> process_cov;
    process_cov <<
        Eigen::Vector<real_t, 3>::Constant(process_noise_trans),
        Eigen::Vector<real_t, 3>::Constant(process_noise_rot),
        Eigen::Vector<real_t, 9>::Constant(process_noise_bias);
    filter.process_noise_covariance = process_cov.asDiagonal();
    std::cout << filter.process_noise_covariance << std::endl;

    filter.measurement_covariance <<
        Eigen::Vector<real_t, 9>::Constant(measurement_noise_hi),
        Eigen::Vector<real_t, 9>::Constant(measurement_noise_lo);
    std::cout << filter.measurement_covariance << std::endl;

    // Iterate
    for (int i = 0; i < coupling_lo.rows(); i++) {
        Low_MeasurementVector meas;
        meas.set_field<High_At_Low_Coupling>(coupling_t(coupling_hi_at_lo(i, Eigen::all)));
        meas.set_field<Low_Coupling>(coupling_t(coupling_lo(i, Eigen::all)));

        // filter.step(dt, meas);
        filter.a_priori_step(dt);
        // std::cout << "1: " << filter.state.get_field<Coupling_Bias>() << std::endl;
        filter.innovation_step(meas);
        // std::cout << "2: " << filter.state.get_field<Coupling_Bias>() << std::endl;
        filter.a_posteriori_step();
        // std::cout << "3: " << filter.state.get_field<Coupling_Bias>() << std::endl;

        biases(i, Eigen::all) = filter.state.get_field<Coupling_Bias>();
        // std::cout << "pose: " << filter.state.get_field<Position>() << filter.state.get_field<Orientation>() << std::endl;
        // std::cout << "bias: " << filter.state.get_field<Coupling_Bias>() << std::endl;
        std::cout << i << std::endl;
    }
    
    return biases;
}