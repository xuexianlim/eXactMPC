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

template <typename T>
T actuatorLen(T q)
{
    T alpha = q(0);
    T beta = q(1);
    T gamma = q(2);

    T lenBoom = 0.00886 * pow(alpha, 4) - 0.0459 * pow(alpha, 3) - 0.0104 * pow(alpha, 2) + 0.2956 * alpha + 1.042;
    T lenArm = 0.0078 * pow(beta, 4) + 0.0917 * pow(beta, 3) + 0.2910 * pow(beta, 2) + 0.0646 * beta + 1.0149;
    T lenBucket = 0.0048 * pow(gamma, 4) + 0.0288 * pow(gamma, 3) + 0.0225 * pow(gamma, 2) - 0.1695 * gamma + 0.9434;

    return vertcat(std::vector<T>{lenBoom, lenArm, lenBucket});
}

template <typename T>
T actuatorVelFactor(T q)
{
    T alpha = q(0);
    T beta = q(1);
    T gamma = q(2);

    T rBoom = 0.0344 * pow(alpha, 3) - 0.1377 * pow(alpha, 2) - 0.0208 * alpha + 0.2956;
    T rArm = 0.0312 * pow(beta, 3) + 0.2751 * pow(beta, 2) + 0.582 * beta + 0.0646;
    T rBucket = 0.0192 * pow(gamma, 3) + 0.0864 * pow(gamma, 2) + 0.045 * gamma - 0.1695;

    return vertcat(std::vector<T>{rBoom, rArm, rBucket});
}

template <typename T>
T actuatorVel(T q, T qDot)
{
    T alpha = q(0);
    T beta = q(1);
    T gamma = q(2);
    T alphaDot = qDot(0);
    T betaDot = qDot(1);
    T gammaDot = qDot(2);

    T lenBoomDot = alphaDot * (0.0344 * pow(alpha, 3) - 0.1377 * pow(alpha, 2) - 0.0208 * alpha + 0.2956);
    T lenArmDot = betaDot * (0.0312 * pow(beta, 3) + 0.2751 * pow(beta, 2) + 0.582 * beta + 0.0646);
    T lenBucketDot = gammaDot * (0.0192 * pow(gamma, 3) + 0.0864 * pow(gamma, 2) + 0.045 * gamma - 0.1695);

    return vertcat(std::vector<T>{lenBoomDot, lenArmDot, lenBucketDot});
}

template <typename T>
T motorVel(T q, T qDot)
{
    T lenDot = actuatorVel(q, qDot);
    T lenBoomDot = lenDot(0);
    T lenArmDot = lenDot(1);
    T lenBucketDot = lenDot(2);

    T angVelMotorBoom = 2444.16 * lenBoomDot;
    T angVelMotorArm = 2444.16 * lenArmDot;
    T angVelMotorBucket = 2444.16 * lenBucketDot;

    return vertcat(std::vector<T>{angVelMotorBoom, angVelMotorArm, angVelMotorBucket});
}

template <typename T>
T motorTorque(T q, T qDot, T qDDot)
{
    T torqueJoints = inverseDynamics(q, qDot, qDDot);
    T torqueBoom = torqueJoints(0);
    T torqueArm = torqueJoints(1);
    T torqueBucket = torqueJoints(2);

    T r = actuatorVelFactor(q);
    T rBoom = r(0);
    T rArm = r(1);
    T rBucket = r(2);

    T torqueMotorBoom = torqueBoom / (2444.16 * rBoom);
    T torqueMotorArm = torqueArm / (2444.16 * rArm);
    T torqueMotorBucket = torqueBucket / (2444.16 * rBucket);

    return vertcat(std::vector<T>{torqueMotorBoom, torqueMotorArm, torqueBucket});
}

#endif