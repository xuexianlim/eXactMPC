import casadi as csd
import excavatorConstants as C

def forwardKinematics(q):
    alpha = q[0]
    beta = q[1]
    gamma = q[2]

    x = C.lenBA*csd.cos(alpha) + C.lenAL*csd.cos(alpha + beta) + C.lenLM*csd.cos(alpha + beta + gamma)
    y = C.lenBA*csd.sin(alpha) + C.lenAL*csd.sin(alpha + beta) + C.lenLM*csd.sin(alpha + beta + gamma)
    theta = alpha + beta + gamma
    return csd.vertcat(x, y, theta)

def inverseKinematics(pose):
    xTip = pose[0]
    yTip = pose[1]
    thetaTip = pose[2]

    xJointBucket = xTip - C.lenLM*csd.cos(thetaTip)
    yJointBucket = yTip - C.lenLM*csd.sin(thetaTip)

    cosBeta = (xJointBucket**2 + yJointBucket**2 - C.lenBA**2 - C.lenAL**2)/(2*C.lenBA*C.lenAL)
    sinBeta = -csd.sqrt(1 - cosBeta**2)

    print(cosBeta, sinBeta)

    sinAlpha = ((C.lenBA + C.lenAL*cosBeta)*yJointBucket - C.lenAL*sinBeta*xJointBucket)/(xJointBucket**2 + yJointBucket**2)
    cosAlpha = ((C.lenBA + C.lenAL*cosBeta)*xJointBucket + C.lenAL*sinBeta*yJointBucket)/(xJointBucket**2 + yJointBucket**2)

    alpha = csd.atan2(sinAlpha, cosAlpha)
    beta = csd.atan2(sinBeta, cosBeta)
    gamma = thetaTip - alpha - beta

    return csd.vertcat(alpha, beta, gamma)

def inverseDynamics(q, qDot, qDDot, F):
    alpha = q[0]
    beta = q[1]
    gamma = q[2]
    alphaDot = qDot[0]
    betaDot = qDot[1]
    gammaDot = qDot[2]

    jacBoom = csd.vertcat(csd.horzcat(-C.lenBCoMBoom*csd.sin(alpha + C.angABCoMBoom), 0, 0),
                          csd.horzcat(C.lenBCoMBoom*csd.cos(alpha + C.angABCoMBoom), 0, 0),
                          csd.horzcat(1, 0, 0))
    jacArm = csd.vertcat(csd.horzcat(-C.lenBA*csd.sin(alpha) - C.lenACoMArm*csd.sin(alpha + beta + C.angLACoMArm),
                                     -C.lenACoMArm*csd.sin(alpha + beta + C.angLACoMArm),
                                     0),
                         csd.horzcat(C.lenBA*csd.cos(alpha) + C.lenACoMArm*csd.cos(alpha + beta + C.angLACoMArm),
                                     C.lenACoMArm*csd.cos(alpha + beta + C.angLACoMArm),
                                     0),
                         csd.horzcat(1, 1, 0))
    jacBucket = csd.vertcat(csd.horzcat(-C.lenBA*csd.sin(alpha) - C.lenAL*csd.sin(alpha + beta) - C.lenLCoMBucket*csd.sin(alpha + beta + gamma + C.angMLCoMBucket),
                                        -C.lenAL*csd.sin(alpha + beta) - C.lenLCoMBucket*csd.sin(alpha + beta + gamma + C.angMLCoMBucket),
                                        -C.lenLCoMBucket*csd.sin(alpha + beta + gamma + C.angMLCoMBucket)),
                            csd.horzcat(C.lenBA*csd.cos(alpha) + C.lenAL*csd.cos(alpha + beta) + C.lenLCoMBucket*csd.cos(alpha + beta + gamma + C.angMLCoMBucket),
                                        C.lenAL*csd.cos(alpha + beta) + C.lenLCoMBucket*csd.cos(alpha + beta + gamma + C.angMLCoMBucket),
                                        C.lenLCoMBucket*csd.cos(alpha + beta + gamma + C.angMLCoMBucket)),
                            csd.horzcat(1, 1, 1))
    jacBoomDot = csd.vertcat(csd.horzcat(-C.lenBCoMBoom*csd.cos(alpha + C.angABCoMBoom)*alphaDot, 0, 0),
                             csd.horzcat(-C.lenBCoMBoom*csd.sin(alpha + C.angABCoMBoom)*alphaDot, 0, 0),
                             csd.horzcat(0, 0, 0))
    jacArmDot = csd.vertcat(csd.horzcat(-C.lenBA*csd.cos(alpha)*alphaDot - C.lenACoMArm*csd.cos(alpha + beta + C.angLACoMArm)*(alphaDot + betaDot),
                                        -C.lenACoMArm*csd.cos(alpha + beta + C.angLACoMArm)*(alphaDot + betaDot),
                                        0),
                            csd.horzcat(-C.lenBA*csd.sin(alpha)*alphaDot - C.lenACoMArm*csd.sin(alpha + beta + C.angLACoMArm)*(alphaDot + betaDot),
                                        -C.lenACoMArm*csd.sin(alpha + beta + C.angLACoMArm)*(alphaDot + betaDot),
                                        0),
                            csd.horzcat(0, 0, 0))
    jacBucketDot = csd.vertcat(csd.horzcat(-C.lenBA*csd.cos(alpha)*alphaDot - C.lenAL*csd.cos(alpha + beta)*(alphaDot + betaDot) - C.lenLCoMBucket*csd.cos(alpha + beta + gamma + C.angMLCoMBucket)*(alphaDot + betaDot + gammaDot),
                                           -C.lenAL*csd.cos(alpha + beta)*(alphaDot + betaDot) - C.lenLCoMBucket*csd.cos(alpha + beta + gamma + C.angMLCoMBucket)*(alphaDot + betaDot + gammaDot),
                                           -C.lenLCoMBucket*csd.cos(alpha + beta + gamma + C.angMLCoMBucket)*(alphaDot + betaDot + gammaDot)),
                               csd.horzcat(-C.lenBA*csd.sin(alpha)*alphaDot - C.lenAL*csd.sin(alpha + beta)*(alphaDot + betaDot) - C.lenLCoMBucket*csd.sin(alpha + beta + gamma + C.angMLCoMBucket)*(alphaDot + betaDot + gammaDot),
                                           -C.lenAL*csd.sin(alpha + beta)*(alphaDot + betaDot) - C.lenLCoMBucket*csd.sin(alpha + beta + gamma + C.angMLCoMBucket)*(alphaDot + betaDot + gammaDot),
                                           -C.lenLCoMBucket*csd.sin(alpha + beta + gamma + C.angMLCoMBucket)*(alphaDot + betaDot + gammaDot)),
                               csd.horzcat(0, 0, 0))
    jacTip = csd.vertcat(csd.horzcat(-C.lenBA*csd.sin(alpha) - C.lenAL*csd.sin(alpha + beta) - C.lenLM*csd.sin(alpha + beta + gamma),
                                     -C.lenAL*csd.sin(alpha + beta) - C.lenLM*csd.sin(alpha + beta + gamma),
                                     -C.lenLM*csd.sin(alpha + beta + gamma)),
                        csd.horzcat(C.lenBA*csd.cos(alpha) + C.lenAL*csd.cos(alpha + beta) + C.lenLM*csd.cos(alpha + beta + gamma),
                                    C.lenAL*csd.cos(alpha + beta) + C.lenLM*csd.cos(alpha + beta + gamma),
                                    C.lenLM*csd.cos(alpha + beta + gamma)),
                        csd.horzcat(1, 1, 1))

    M = (csd.transpose(jacBoom[0:2, :])@C.massBoom@jacBoom[0:2, :] + csd.transpose(jacBoom[2, :])@C.moiBoom@jacBoom[2, :] +
         csd.transpose(jacArm[0:2, :])@C.massArm@jacArm[0:2, :] + csd.transpose(jacArm[2, :])@C.moiArm@jacArm[2, :] +
         csd.transpose(jacBucket[0:2, :])@C.massBucket@jacBucket[0:2, :] + csd.transpose(jacBucket[2, :])@C.moiBucket@jacBucket[2, :])
    b = (csd.transpose(jacBoom[0:2, :])@C.massBoom@jacBoomDot[0:2, :]@qDot + csd.transpose(jacBoom[2, :])@(C.moiBoom@jacBoomDot[2, :]@qDot) +
         csd.transpose(jacArm[0:2, :])@C.massArm@jacArmDot[0:2, :]@qDot + csd.transpose(jacArm[2, :])@(C.moiArm@jacArmDot[2, :]@qDot) +
         csd.transpose(jacBucket[0:2, :])@C.massBucket@jacBucketDot[0:2, :]@qDot + csd.transpose(jacBucket[2, :])@(C.moiBucket@jacBucketDot[2, :]@qDot))
    g = (-csd.transpose(jacBoom[0:2, :])@C.massBoom@C.g -
         csd.transpose(jacArm[0:2, :])@C.massArm@C.g -
         csd.transpose(jacBucket[0:2, :])@C.massBucket@C.g)
    extF = csd.transpose(jacTip[0:2, :])@F

    return M@qDDot + b + g - extF

def actuatorLen(q):
    alpha = q[0]
    beta = q[1]
    gamma = q[2]

    # R = 2*C.lenBC
    # theta = csd.atan2(C.iBC[1], C.iBC[0])
    # lenBoom = csd.sqrt(-R*C.lenBD*csd.cos(alpha - theta + C.angABD) + C.lenBD**2 + C.lenBC**2)
    lenBoom = 0.0086*alpha**4 - 0.0459*alpha**3 - 0.0104*alpha**2 + 0.2956*alpha + 1.042

    # R = 2*C.lenAE
    # theta = csd.atan2(C.bAE[0], C.bAE[1])
    # lenArm = csd.sqrt(-R*C.lenAF*csd.sin(beta + theta + C.angFAL) + C.lenAF**2 + C.lenAE**2)
    lenArm = 0.0078*beta**4 + 0.0917*beta**3 + 0.2910*beta**2 + 0.0646*beta + 1.0149

    lenBucket = 0.0048*gamma**4 + 0.0288*gamma**3 + 0.0225*gamma**2 - 0.1695*gamma + 0.9434

    return csd.vertcat(lenBoom, lenArm, lenBucket)

def actuatorVelFactor(q):
    alpha = q[0]
    beta = q[1]
    gamma = q[2]

    velFactorBoom = 0.0344*alpha**3 - 0.1377*alpha**2 - 0.0208*alpha + 0.2956
    velFactorArm = 0.0312*beta**3 + 0.2751*beta**2 + 0.582*beta + 0.0646
    velFactorBucket = 0.0192*gamma**3 + 0.0864*gamma**2 + 0.045*gamma - 0.1695

    return csd.vertcat(velFactorBoom, velFactorArm, velFactorBucket)

def actuatorVel(q, qDot):
    alpha = q[0]
    beta = q[1]
    gamma = q[2]
    alphaDot = qDot[0]
    betaDot = qDot[1]
    gammaDot = qDot[2]

    # len = actuatorLen(q)
    # lenBoom = len[0]
    # lenArm = len[1]
    # lenBucket = len[2]

    # angBDC = csd.acos((lenBoom**2 + C.lenBD**2 - C.lenBC**2)/(2*lenBoom*C.lenBD))
    # theta = csd.pi/2 - angBDC
    # lenBoomDot = alphaDot*C.lenBD*csd.cos(theta)
    lenBoomDot = alphaDot*(0.0344*alpha**3 - 0.1377*alpha**2 - 0.0208*alpha + 0.2956)

    # angAFE = csd.acos((lenArm**2 + C.lenAF**2 - C.lenAE**2)/(2*lenArm*C.lenAF))
    # theta = angAFE - csd.pi/2
    # lenArmDot = -betaDot*C.lenAF*csd.cos(theta)
    lenArmDot = betaDot*(0.0312*beta**3 + 0.2751*beta**2 + 0.582*beta + 0.0646)

    lenBucketDot = gammaDot*(0.0192*gamma**3 + 0.0864*gamma**2 + 0.045*gamma - 0.1695)

    return csd.vertcat(lenBoomDot, lenArmDot, lenBucketDot)

def motorVel(q, qDot):
    lenDot = actuatorVel(q, qDot)
    lenBoomDot = lenDot[0]
    lenArmDot = lenDot[1]
    lenBucketDot = lenDot[2]

    angVelMotorBoom = 2444.16*lenBoomDot
    angVelMotorArm = 2444.16*lenArmDot
    angVelMotorBucket = 2444.16*lenBucketDot

    return csd.vertcat(angVelMotorBoom, angVelMotorArm, angVelMotorBucket)

def motorTorque(q, qDot, qDDot, F):
    T = inverseDynamics(q, qDot, qDDot, F)
    TBoom = T[0]
    TArm = T[1]
    TBucket = T[2]

    r = actuatorVelFactor(q)
    rBoom = r[0]
    rArm = r[1]
    rBucket = r[2]

    TMotorBoom = TBoom/(2444.16*rBoom)
    TMotorArm = TArm/(2444.16*rArm)
    TMotorBucket = TBucket/(2444.16*rBucket)

    return csd.vertcat(TMotorBoom, TMotorArm, TMotorBucket)