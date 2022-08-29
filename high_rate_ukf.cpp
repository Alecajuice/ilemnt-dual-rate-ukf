#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#define UKF_DOUBLE_PRECISION
// #define SQRT_CORE
// #define SYSTEM_DYNAMICS
#include "UKF/Types.h"
#include "UKF/Integrator.h"
#include "UKF/StateVector.h"
#include "UKF/MeasurementVector.h"
#include "UKF/Core.h"
#include "dual_rate_ukf.h"
#include <iostream>
#include <chrono>

/* State vector definition */
enum High_States {
    Position,
#ifdef SYSTEM_DYNAMICS
    Velocity,
    Acceleration,
    Angular_Velocity,
#endif
    Orientation
};

using High_StateVector = UKF::StateVector<
    UKF::Field<Position, UKF::Vector<3>>,
#ifdef SYSTEM_DYNAMICS
    UKF::Field<Velocity, UKF::Vector<3>>,
    UKF::Field<Acceleration, UKF::Vector<3>>,
    UKF::Field<Angular_Velocity, UKF::Vector<3>>,
#endif
    UKF::Field<Orientation, UKF::Vector<3>>
>;

/* State transition function */
namespace UKF {

/* SFWA state vector process model. */
template <> template <>
High_StateVector High_StateVector::derivative<>() const {
    High_StateVector output;
    output.set_field<Position>(UKF::Vector<3>(0, 0, 0));
    output.set_field<Orientation>(UKF::Vector<3>(0, 0, 0));
#ifdef SYSTEM_DYNAMICS
    output.set_field<Position>(get_field<Velocity>());
    output.set_field<Velocity>(get_field<Acceleration>());
    output.set_field<Acceleration>(UKF::Vector<3>(0, 0, 0));
    output.set_field<Orientation>(get_field<Angular_Velocity>());
    output.set_field<Angular_Velocity>(UKF::Vector<3>(0, 0, 0));
#endif
    return output;
}

}

/* Measurement vector definition */
enum High_Measurements {
    High_Coupling
};

using High_MeasurementVector = UKF::DynamicMeasurementVector<
    UKF::Field<High_Coupling, coupling_t>
>;

/* Measurement functions */
namespace UKF {

template <> template <>
coupling_t High_MeasurementVector::expected_measurement
<High_StateVector, High_Coupling>(const High_StateVector& state) {
    auto fk = forward_kinematics(state.get_field<Position>(),
                                 state.get_field<Orientation>(),
                                 calibration_hi);
    return fk;
}

}

/* UKF Core definition */
#ifndef SQRT_CORE
using High_Core = UKF::Core<
#else
using High_Core = UKF::SquareRootParameterEstimationCore<
#endif
    High_StateVector,
    High_MeasurementVector,
    UKF::IntegratorRK4
>;

/* Run UKF */
Eigen::Matrix<real_t, Eigen::Dynamic, 6> run_high_ukf(
    Eigen::Matrix<real_t, Eigen::Dynamic, 9> &coupling_hi)
{
    // Return value
    Eigen::Matrix<real_t, Eigen::Dynamic, 6> poses(coupling_hi.rows(), 6);

    // Const values
#ifndef SQRT_CORE
    const real_t initial_variance = 1e-6;
#else
    const real_t initial_variance = 1e-3;
#endif

    const real_t base_process_noise = 1e-5;
    const real_t process_noise_trans = base_process_noise;
    const real_t process_noise_rot = base_process_noise;
    const real_t process_noise_vel = base_process_noise * 10;
    const real_t process_noise_acc = base_process_noise * 100;
    const real_t process_noise_w = base_process_noise * 10;

    const real_t measurement_noise = 1e-7;

    const real_t dt = 1.0/1500;

    // Initialize filter
    High_Core filter;

    // TODO: either solve for initial pose or get it from matlab
    // filter.state.set_field<Position>(UKF::Vector<3>({0.1107, 0.1600, 0.1177})); // noise_still_30s
    // filter.state.set_field<Orientation>(UKF::Vector<3>({-2.2159, 0.8268, -0.3796})); // noise_still_30s
    filter.state.set_field<Position>(UKF::Vector<3>({0.0883, 0.1543, 0.0977})); // aluminum_sheet
    filter.state.set_field<Orientation>(UKF::Vector<3>({-0.9552, 2.1783, 0.4121})); // aluminum_sheet
#ifdef SYSTEM_DYNAMICS
    filter.state.set_field<Velocity>(UKF::Vector<3>(0, 0, 0));
    filter.state.set_field<Acceleration>(UKF::Vector<3>(0, 0, 0));
    filter.state.set_field<Angular_Velocity>(UKF::Vector<3>(0, 0, 0));
#endif

#ifndef SQRT_CORE
    filter.covariance = High_StateVector::CovarianceMatrix::Identity() * initial_variance;
#else
    filter.root_covariance = High_StateVector::CovarianceMatrix::Identity() * initial_variance;
#endif

    Eigen::Vector<real_t, 6> process_cov;
    process_cov <<
        Eigen::Vector<real_t, 3>::Constant(process_noise_trans),
#ifdef SYSTEM_DYNAMICS
        Eigen::Vector<real_t, 3>::Constant(process_noise_vel),
        Eigen::Vector<real_t, 3>::Constant(process_noise_acc),
        Eigen::Vector<real_t, 3>::Constant(process_noise_w),
#endif
        Eigen::Vector<real_t, 3>::Constant(process_noise_rot);
#ifndef SQRT_CORE
    process_cov = process_cov.array().pow(2);
    filter.process_noise_covariance = process_cov.asDiagonal();
#else
    filter.process_noise_root_covariance = process_cov.asDiagonal();
#endif

    Eigen::Vector<real_t, 9> measurement_cov;
    measurement_cov <<
        Eigen::Vector<real_t, 9>::Constant(measurement_noise);
#ifndef SQRT_CORE
    measurement_cov = measurement_cov.array().pow(2);
    filter.measurement_covariance = measurement_cov;
#else
    filter.measurement_root_covariance = measurement_cov;
#endif

    // Iterate
    auto starttime = std::chrono::high_resolution_clock::now();
    auto prevtime = starttime;
    int num_iters = coupling_hi.rows();
    for (int i = 0; i < num_iters; i++) {
        std::cout << "iteration: " << i << std::endl;

        High_MeasurementVector meas;
        coupling_t hi = coupling_t(coupling_hi(i, Eigen::all));
        meas.set_field<High_Coupling>(hi);

        // filter.step(dt, meas);
#ifndef SQRT_CORE
        filter.a_priori_step(dt);
#else
        filter.a_priori_step();
#endif
        filter.innovation_step(meas);
        filter.a_posteriori_step();

        poses(i, Eigen::seq(0, 2)) = filter.state.get_field<Position>();
        poses(i, Eigen::seq(3, 5)) = filter.state.get_field<Orientation>();
        auto curtime = std::chrono::high_resolution_clock::now();
        // std::cout << "total time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(curtime - prevtime).count() << std::endl;
        prevtime = curtime;
    }

    auto endtime = std::chrono::high_resolution_clock::now();
    auto avgtime = std::chrono::duration_cast<std::chrono::nanoseconds>((endtime - starttime) / num_iters).count();
    std::cout << "average time per iteration: " << avgtime << std::endl;
    
    return poses;
}