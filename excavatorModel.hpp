#ifndef EXCAVATOR_MODEL_H
#define EXCAVATOR_MODEL_H

#include "casadi/casadi.hpp"

using namespace casadi;

DM forwardKinematics(DM q);
DM inverseDynamics(DM q, DM qDot, DM qDDot);
DM angle2Length(DM q);
DM length2Angle(DM l);
DM jointTorque2ActuatorForce(DM T);
DM actuatorForce2MotorTorque(DM F);

#endif