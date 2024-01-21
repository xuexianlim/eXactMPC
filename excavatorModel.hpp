#ifndef EXCAVATOR_MODEL_H
#define EXCAVATOR_MODEL_H

#include <vector>
#include "casadi/casadi.hpp"
#include "excavatorConstants.hpp"

using namespace casadi;

template <typename T>
T forwardKinematics(T q)
{
    T alpha = q(0);
    T beta = q(1);
    T gamma = q(2);

    T x = lenBA * cos(alpha) + lenAL * cos(alpha + beta) + lenLM * cos(alpha + beta + gamma);
    T y = lenBA * sin(alpha) + lenAL * sin(alpha + beta) + lenLM * sin(alpha + beta + gamma);
    T theta = alpha + beta + gamma;

    return vertcat(std::vector<T>{x, y, theta});
}

template <typename T>
T inverseDynamics(T q, T qDot, T qDDot)
{
    T alpha = q(0);
    T beta = q(1);
    T gamma = q(2);
    T alphaDot = qDot(0);
    T betaDot = qDot(1);
    T gammaDot = qDot(2);

    T jacBoomPos = vertcat(std::vector<T>{horzcat(std::vector<T>{-lenBCoMBoom * sin(alpha + angABCoMBoom), 0, 0}),
                                          horzcat(std::vector<T>{lenBCoMBoom * cos(alpha + angABCoMBoom), 0, 0})});
    T jacBoomRot = horzcat(std::vector<T>{1, 0, 0});

    T jacArmPos = vertcat(std::vector<T>{horzcat(std::vector<T>{-lenBA * sin(alpha) - lenACoMArm * sin(alpha + beta + angLACoMArm),
                                                                -lenACoMArm * sin(alpha + beta + angLACoMArm),
                                                                0}),
                                         horzcat(std::vector<T>{lenBA * cos(alpha) + lenACoMArm * cos(alpha + beta + angLACoMArm),
                                                                lenACoMArm * cos(alpha + beta + angLACoMArm),
                                                                0})});
    T jacArmRot = horzcat(std::vector<T>{1, 1, 0});

    T jacBucketPos = vertcat(std::vector<T>{horzcat(std::vector<T>{-lenBA * sin(alpha) - lenAL * sin(alpha + beta) - lenLCoMBucket * sin(alpha + beta + gamma + angMLCoMBucket),
                                                                   -lenAL * sin(alpha + beta) - lenLCoMBucket * sin(alpha + beta + gamma + angMLCoMBucket),
                                                                   -lenLCoMBucket * sin(alpha + beta + gamma + angMLCoMBucket)}),
                                            horzcat(std::vector<T>{lenBA * cos(alpha) + lenAL * cos(alpha + beta) + lenLCoMBucket * cos(alpha + beta + gamma + angMLCoMBucket),
                                                                   lenAL * cos(alpha + beta) + lenLCoMBucket * cos(alpha + beta + gamma + angMLCoMBucket),
                                                                   lenLCoMBucket * cos(alpha + beta + gamma + angMLCoMBucket)})});
    T jacBucketRot = horzcat(std::vector<T>{1, 1, 1});

    T jacBoomDotPos = vertcat(std::vector<T>{horzcat(std::vector<T>{-lenBCoMBoom * cos(alpha + angABCoMBoom) * alphaDot, 0, 0}),
                                             horzcat(std::vector<T>{-lenBCoMBoom * sin(alpha + angABCoMBoom) * alphaDot, 0, 0})});
    T jacBoomDotRot = horzcat(std::vector<T>{0, 0, 0});

    T jacArmDotPos = vertcat(std::vector<T>{horzcat(std::vector<T>{-lenBA * cos(alpha) * alphaDot - lenACoMArm * cos(alpha + beta + angLACoMArm) * (alphaDot + betaDot),
                                                                   -lenACoMArm * cos(alpha + beta + angLACoMArm) * (alphaDot + betaDot),
                                                                   0}),
                                            horzcat(std::vector<T>{-lenBA * sin(alpha) * alphaDot - lenACoMArm * sin(alpha + beta + angLACoMArm) * (alphaDot + betaDot),
                                                                   -lenACoMArm * sin(alpha + beta + angLACoMArm) * (alphaDot + betaDot),
                                                                   0})});
    T jacArmDotRot = horzcat(std::vector<T>{0, 0, 0});

    T jacBucketDotPos = vertcat(std::vector<T>{horzcat(std::vector<T>{-lenBA * cos(alpha) * alphaDot - lenAL * cos(alpha + beta) * (alphaDot + betaDot) - lenLCoMBucket * cos(alpha + beta + gamma + angMLCoMBucket) * (alphaDot + betaDot + gammaDot),
                                                                      -lenAL * cos(alpha + beta) * (alphaDot + betaDot) - lenLCoMBucket * cos(alpha + beta + gamma + angMLCoMBucket) * (alphaDot + betaDot + gammaDot),
                                                                      -lenLCoMBucket * cos(alpha + beta + gamma + angMLCoMBucket) * (alphaDot + betaDot + gammaDot)}),
                                               horzcat(std::vector<T>{-lenBA * sin(alpha) * alphaDot - lenAL * sin(alpha + beta) * (alphaDot + betaDot) - lenLCoMBucket * sin(alpha + beta + gamma + angMLCoMBucket) * (alphaDot + betaDot + gammaDot),
                                                                      -lenAL * sin(alpha + beta) * (alphaDot + betaDot) - lenLCoMBucket * sin(alpha + beta + gamma + angMLCoMBucket) * (alphaDot + betaDot + gammaDot),
                                                                      -lenLCoMBucket * sin(alpha + beta + gamma + angMLCoMBucket) * (alphaDot + betaDot + gammaDot)})});
    T jacBucketDotRot = horzcat(std::vector<T>{0, 0, 0});

    T M = mtimes(jacBoomPos.T() * massBoom, jacBoomPos) + mtimes(jacBoomRot.T() * moiBoom, jacBoomRot) +
          mtimes(jacArmPos.T() * massArm, jacArmPos) + mtimes(jacArmRot.T() * moiArm, jacArmRot) +
          mtimes(jacBucketPos.T() * massBucket, jacBucketPos) + mtimes(jacBucketRot.T() * moiBucket, jacBucketRot);
    T b = mtimes(jacBoomPos.T() * massBoom, mtimes(jacBoomDotPos, qDot)) + mtimes(jacBoomRot.T() * moiBoom, mtimes(jacBoomDotRot, qDot)) +
          mtimes(jacArmPos.T() * massArm, mtimes(jacArmDotPos, qDot)) + mtimes(jacArmRot.T() * moiArm, mtimes(jacArmDotRot, qDot)) +
          mtimes(jacBucketPos.T() * massBucket, mtimes(jacBucketDotPos, qDot)) + mtimes(jacBucketRot.T() * moiBucket, mtimes(jacBucketDotRot, qDot));
    T g = -mtimes(jacBoomPos.T() * massBoom, gravity) - mtimes(jacArmPos.T() * massArm, gravity) - mtimes(jacBucketPos.T() * massBucket, gravity);

    return mtimes(M, qDDot) + b + g;
}

DM angle2Length(DM q);
DM length2Angle(DM l);
DM jointTorque2ActuatorForce(DM T);
DM actuatorForce2MotorTorque(DM F);

#endif