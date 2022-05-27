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