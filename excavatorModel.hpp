#ifndef EXCAVATOR_MODEL_H
#define EXCAVATOR_MODEL_H

#include <vector>
#include "casadi/casadi.hpp"
#include "excavatorConstants.hpp"

using namespace casadi;

enum Mode
{
    NO_LOAD,
    LIFT,
    DIG
};

enum DutyCycle
{
    S1,
    S2_60,
    S2_30,
    PEAK
};

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
T inverseKinematics(T pose)
{
    T xTip = pose(0);
    T yTip = pose(1);
    T thetaTip = pose(2);

    T xJointBucket = xTip - lenLM * cos(thetaTip);
    T yJointBucket = yTip - lenLM * sin(thetaTip);

    T cosBeta = (pow(xJointBucket, 2) + pow(yJointBucket, 2) - pow(lenBA, 2) - pow(lenAL, 2)) / (2 * lenBA * lenAL);
    T sinBeta = -sqrt(1 - pow(cosBeta, 2));

    T sinAlpha = ((lenBA + lenAL * cosBeta) * yJointBucket - lenAL * sinBeta * xJointBucket) / (pow(xJointBucket, 2) + pow(yJointBucket, 2));
    T cosAlpha = ((lenBA + lenAL * cosBeta) * xJointBucket + lenAL * sinBeta * yJointBucket) / (pow(xJointBucket, 2) + pow(yJointBucket, 2));

    T alpha = atan2(sinAlpha, cosAlpha);
    T beta = atan2(sinBeta, cosBeta);
    T gamma = thetaTip - alpha - beta;

    return vertcat(std::vector<T>{alpha, beta, gamma});
}

template <typename T>
T inverseDynamics(T q, T qDot, T qDDot, T F)
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
    T extF = mtimes(jacBucketPos.T(), F);

    return mtimes(M, qDDot) + b + g - extF;
}

template <typename T>
T jointAngles(T len)
{
    T lenBoom = len(0);
    T lenArm = len(1);
    T lenBucket = len(2);

    T R = lenBC;
    T theta = atan2(iBC(1), iBC(0));
    T alpha = acos((-pow(lenBoom, 2) + pow(lenBD, 2) + pow(lenBC, 2)) / (2 * R * lenBD)) + theta - angABD;

    R = lenAE;
    theta = atan2(bAE(0), bAE(1));
    T beta = asin((-pow(lenArm, 2) + pow(lenAF, 2) + pow(lenAE, 2)) / (2 * R * lenAF)) - theta - angFAL;

    R = lenJG;
    theta = atan2(aJG(0), aJG(1));
    T angLJH = asin((-pow(lenBucket, 2) + pow(lenHJ, 2) + pow(lenJG, 2)) / (2 * R * lenHJ)) - theta;

    R = sqrt(pow((lenJL - lenHJ * cos(angLJH)), 2) + pow((lenHJ * sin(angLJH)), 2));
    theta = atan2(lenHJ * sin(angLJH), lenJL - lenHJ * cos(angLJH));
    T gamma = acos((pow(lenHK, 2) - pow(lenJL, 2) - pow(lenLK, 2) - pow(lenHJ, 2) + 2 * lenJL * lenHJ * cos(angLJH)) / (2 * R * lenLK)) - theta - angKLM;

    return vertcat(std::vector<T>{alpha, beta, gamma});
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
T motorTorque(T q, T qDot, T qDDot, T F)
{
    T torqueJoints = inverseDynamics(q, qDot, qDDot, F);
    T torqueBoom = torqueJoints(0);
    T torqueArm = torqueJoints(1);
    T torqueBucket = torqueJoints(2);

    T r = actuatorVelFactor(q);
    T rBoom = r(0);
    T rArm = r(1);
    T rBucket = r(2);

    T torqueMotorBoom = torqueBoom / (2444.16 * rBoom * 0.64);
    T torqueMotorArm = torqueArm / (2444.16 * rArm * 0.64);
    T torqueMotorBucket = torqueBucket / (2444.16 * rBucket * 0.64);

    return vertcat(std::vector<T>{torqueMotorBoom, torqueMotorArm, torqueMotorBucket});
}

template <typename T>
T motorTorqueLimit(T motorVel, DutyCycle dutyCycle)
{
    switch (dutyCycle)
    {
    case S1:
        return -1.4073e-7 * pow(motorVel, 3) + 1.7961e-5 * pow(motorVel, 2) - 0.0147 * motorVel + 19.9091;
    case S2_60:
        return -1.3568e-10 * pow(motorVel, 4) - 1.4682e-7 * pow(motorVel, 3) + 5.5744e-5 * pow(motorVel, 2) - 0.0159 * motorVel + 25.0769;
    case S2_30:
        return -2.0433e-7 * pow(motorVel, 3) + 5.1865e-5 * pow(motorVel, 2) - 0.0105 * motorVel + 30.9119;
    case PEAK:
        return 4.0269e-9 * pow(motorVel, 4) - 3.7090e-6 * pow(motorVel, 3) + 8.513e-4 * pow(motorVel, 2) - 0.05787 * motorVel + 60.2614;
    default:
        return -1.4073e-7 * pow(motorVel, 3) + 1.7961e-5 * pow(motorVel, 2) - 0.0147 * motorVel + 19.9091;
    }
}

#endif