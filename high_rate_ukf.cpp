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
#include <chrono>

/* State vector definition */
enum High_States {
    Position,
    // Velocity,
    // Acceleration,
    Orientation,
    // Angular_Velocity
};

using High_StateVector = UKF::StateVector<
    UKF::Field<Position, UKF::Vector<3>>,
    // UKF::Field<Velocity, UKF::Vector<3>>,
    // UKF::Field<Acceleration, UKF::Vector<3>>,
    UKF::Field<Orientation, UKF::Vector<3>>
    // UKF::Field<Angular_Velocity, UKF::Vector<3>>
>;

/* State transition function */
namespace UKF {

/* SFWA state vector process model. */
template <> template <>
High_StateVector High_StateVector::derivative<>() const {
    High_StateVector output;
    output.set_field<Position>(UKF::Vector<3>(0, 0, 0));
    output.set_field<Orientation>(UKF::Vector<3>(0, 0, 0));
    // output.set_field<Position>(get_field<Velocity>());
    // output.set_field<Velocity>(get_field<Acceleration>());
    // output.set_field<Acceleration>(UKF::Vector<3>(0, 0, 0));
    // output.set_field<Orientation>(get_field<Angular_Velocity>());
    // output.set_field<Angular_Velocity>(UKF::Vector<3>(0, 0, 0));
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
using High_Core = UKF::Core< // Non-Sqrt Core
// using High_Core = UKF::SquareRootCore< // Sqrt Core
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
    const real_t initial_variance = 1e-6; // Non-Sqrt Core
    // const real_t initial_variance = 1e-3; // Sqrt Core

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

    // filter.state.set_field<Position>(UKF::Vector<3>({0.1107, 0.1600, 0.1177})); // noise_still_30s
    // filter.state.set_field<Orientation>(UKF::Vector<3>({-2.2159, 0.8268, -0.3796})); // noise_still_30s
    filter.state.set_field<Position>(UKF::Vector<3>({0.0883, 0.1543, 0.0977})); // aluminum_sheet
    filter.state.set_field<Orientation>(UKF::Vector<3>({-0.9552, 2.1783, 0.4121})); // aluminum_sheet
    // filter.state.set_field<Velocity>(UKF::Vector<3>(0, 0, 0));
    // filter.state.set_field<Acceleration>(UKF::Vector<3>(0, 0, 0));
    // filter.state.set_field<Angular_Velocity>(UKF::Vector<3>(0, 0, 0));

    filter.covariance = High_StateVector::CovarianceMatrix::Identity() * initial_variance; // Non-Sqrt Core
    // filter.root_covariance = High_StateVector::CovarianceMatrix::Identity() * initial_variance; // Sqrt Core

    Eigen::Vector<real_t, 6> process_cov;
    process_cov <<
        Eigen::Vector<real_t, 3>::Constant(process_noise_trans),
        // Eigen::Vector<real_t, 3>::Constant(process_noise_vel),
        // Eigen::Vector<real_t, 3>::Constant(process_noise_acc),
        Eigen::Vector<real_t, 3>::Constant(process_noise_rot);
        // Eigen::Vector<real_t, 3>::Constant(process_noise_w);
    process_cov = process_cov.array().pow(2); // Non-Sqrt Core
    filter.process_noise_covariance = process_cov.asDiagonal(); // Non-Sqrt Core
    // filter.process_noise_root_covariance = process_cov.asDiagonal(); // Sqrt Core

    Eigen::Vector<real_t, 9> measurement_cov;
    measurement_cov <<
        Eigen::Vector<real_t, 9>::Constant(measurement_noise);
    measurement_cov = measurement_cov.array().pow(2); // Non-Sqrt Core
    filter.measurement_covariance = measurement_cov; // Non-Sqrt Core
    // filter.measurement_root_covariance = measurement_cov; // Sqrt Core

    // Iterate
    auto prevtime = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < 1000; i++) {
    // for (int i = 0; i < coupling_hi.rows(); i++) {
        std::cout << "iteration: " << i << std::endl;

        High_MeasurementVector meas;
        coupling_t hi = coupling_t(coupling_hi(i, Eigen::all));
        meas.set_field<High_Coupling>(hi);

        // filter.step(dt, meas);
        filter.a_priori_step(dt);
        auto time1 = std::chrono::high_resolution_clock::now();
        // std::cout << "1: " << filter.state.get_field<Coupling_Bias>() << std::endl;
        filter.innovation_step(meas);
        auto time2 = std::chrono::high_resolution_clock::now();
        // std::cout << "2: " << filter.state.get_field<Coupling_Bias>() << std::endl;
        filter.a_posteriori_step();
        auto time3 = std::chrono::high_resolution_clock::now();
        // std::cout << "3: " << filter.state.get_field<Coupling_Bias>() << std::endl;

        poses(i, Eigen::seq(0, 2)) = filter.state.get_field<Position>();
        poses(i, Eigen::seq(3, 5)) = filter.state.get_field<Orientation>();
        // std::cout << "time1: " << std::chrono::duration_cast<std::chrono::nanoseconds>(time1 - prevtime).count() << std::endl;
        // std::cout << "time2: " << std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time1).count() << std::endl;
        // std::cout << "time3: " << std::chrono::duration_cast<std::chrono::nanoseconds>(time3 - time2).count() << std::endl;
        auto curtime = std::chrono::high_resolution_clock::now();
        std::cout << "total time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(curtime - prevtime).count() << std::endl;
        prevtime = curtime;
    }

    // std::cout << poses << std::endl;
    
    return poses;
}