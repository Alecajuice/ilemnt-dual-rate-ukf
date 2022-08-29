#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#define UKF_DOUBLE_PRECISION
// #define SQRT_CORE
#include "UKF/Types.h"
#include "UKF/Integrator.h"
#include "UKF/StateVector.h"
#include "UKF/MeasurementVector.h"
#include "UKF/Core.h"
#include "dual_rate_ukf.h"
#include <iostream>
#include <chrono>

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
    auto fk = forward_kinematics(state.get_field<Position>(),
                                 state.get_field<Orientation>(),
                                 calibration_hi);
    auto bias = state.get_field<Coupling_Bias>();
    auto retval = fk + bias;
    return retval;
}

template <> template <>
coupling_t Low_MeasurementVector::expected_measurement
<Low_StateVector, Low_Coupling>(const Low_StateVector& state) {
    auto retval = forward_kinematics(state.get_field<Position>(),
                              state.get_field<Orientation>(),
                              calibration_lo);
    return retval;
}

}

/* UKF Core definition */
#ifndef SQRT_CORE
using Low_Core = UKF::Core<
#else
using Low_Core = UKF::SquareRootCore<
#endif
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
    Eigen::Matrix<real_t, Eigen::Dynamic, 3> positions(coupling_lo.rows(), 3);
    Eigen::Matrix<real_t, Eigen::Dynamic, 9> biases(coupling_lo.rows(), 9);

    // Const values
#ifndef SQRT_CORE
    const real_t initial_variance = 1e-6;
#else
    const real_t initial_variance = 1e-3;
#endif

    const real_t base_process_noise = 2e-3;
    const real_t process_noise_trans = base_process_noise;
    const real_t process_noise_rot = base_process_noise * 2;
    const real_t process_noise_bias = 1e-5;

    const real_t measurement_noise_lo = 1e-7;
    const real_t measurement_noise_hi = measurement_noise_lo * 100;

    const real_t dt = 1.0/1500; // TODO: dt is wrong but actually it doesn't matter

    // Initialize filter
    Low_Core filter;

    // filter.state.set_field<Position>(UKF::Vector<3>({0.1107, 0.1600, 0.1177})); // noise_still_30s
    // filter.state.set_field<Orientation>(UKF::Vector<3>({-2.2159, 0.8268, -0.3796})); // noise_still_30s
    filter.state.set_field<Position>(UKF::Vector<3>({0.0883, 0.1543, 0.0977})); // aluminum_sheet
    filter.state.set_field<Orientation>(UKF::Vector<3>({-0.9552, 2.1783, 0.4121})); // aluminum_sheet
    filter.state.set_field<Coupling_Bias>(UKF::Vector<9>::Zero());

#ifndef SQRT_CORE
    filter.covariance = Low_StateVector::CovarianceMatrix::Identity() * initial_variance;
#else
    filter.root_covariance = Low_StateVector::CovarianceMatrix::Identity() * initial_variance;
#endif

    Eigen::Vector<real_t, 15> process_cov;
    process_cov <<
        Eigen::Vector<real_t, 3>::Constant(process_noise_trans),
        Eigen::Vector<real_t, 3>::Constant(process_noise_rot),
        Eigen::Vector<real_t, 9>::Constant(process_noise_bias);
#ifndef SQRT_CORE
    process_cov = process_cov.array().pow(2);
    filter.process_noise_covariance = process_cov.asDiagonal();
#else
    filter.process_noise_root_covariance = process_cov.asDiagonal();
#endif

    Eigen::Vector<real_t, 18> measurement_cov;
    measurement_cov <<
        Eigen::Vector<real_t, 9>::Constant(measurement_noise_hi),
        Eigen::Vector<real_t, 9>::Constant(measurement_noise_lo);
#ifndef SQRT_CORE
    measurement_cov = measurement_cov.array().pow(2);
    filter.measurement_covariance = measurement_cov;
#else
    filter.measurement_root_covariance = measurement_cov;
#endif

    // Iterate
    auto prevtime = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < coupling_lo.rows(); i++) {
        std::cout << "iteration: " << i << std::endl;

        Low_MeasurementVector meas;
        coupling_t hi_at_lo = coupling_t(coupling_hi_at_lo(i, Eigen::all));
        coupling_t lo = coupling_t(coupling_lo(i, Eigen::all));
        meas.set_field<High_At_Low_Coupling>(hi_at_lo);
        meas.set_field<Low_Coupling>(lo);

        // filter.step(dt, meas);
        filter.a_priori_step(dt);
        filter.innovation_step(meas);
        filter.a_posteriori_step();

        positions(i, Eigen::all) = filter.state.get_field<Position>();
        biases(i, Eigen::all) = filter.state.get_field<Coupling_Bias>();
        auto curtime = std::chrono::high_resolution_clock::now();
        std::cout << "total time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(curtime - prevtime).count() << std::endl;
        prevtime = curtime;
    }
    
    return biases;
}