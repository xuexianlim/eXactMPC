#include <vector>
#include "casadi/casadi.hpp"
#include "excavatorConstants.hpp"
#include "utils.hpp"

using namespace casadi;

DM forwardKinematics(DM q)
{
    DM alpha = q(0);
    DM beta = q(1);
    DM gamma = q(2);

    DM x = lenBA * cos(alpha) + lenAL * cos(alpha + beta) + lenLM * cos(alpha + beta + gamma);
    DM y = lenBA * sin(alpha) + lenAL * sin(alpha + beta) + lenLM * sin(alpha + beta + gamma);
    DM theta = alpha + beta + gamma;

    return vertcat(std::vector<DM>{x, y, theta});
}

DM inverseDynamics(DM q, DM qDot, DM qDDot)
{
    DM alpha = q(0);
    DM beta = q(1);
    DM gamma = q(2);
    DM alphaDot = qDot(0);
    DM betaDot = qDot(1);
    DM gammaDot = qDot(2);

    DM jacBoomPos = vertcat(std::vector<DM>{horzcat(std::vector<DM>{-lenBCoMBoom * sin(alpha + angABCoMBoom), 0, 0}),
                                            horzcat(std::vector<DM>{lenBCoMBoom * cos(alpha + angABCoMBoom), 0, 0})});
    DM jacBoomRot = horzcat(std::vector<DM>{1, 0, 0});

    DM jacArmPos = vertcat(std::vector<DM>{horzcat(std::vector<DM>{-lenBA * sin(alpha) - lenACoMArm * sin(alpha + beta + angLACoMArm),
                                                                   -lenACoMArm * sin(alpha + beta + angLACoMArm),
                                                                   0}),
                                           horzcat(std::vector<DM>{lenBA * cos(alpha) + lenACoMArm * cos(alpha + beta + angLACoMArm),
                                                                   lenACoMArm * cos(alpha + beta + angLACoMArm),
                                                                   0})});
    DM jacArmRot = horzcat(std::vector<DM>{1, 1, 0});

    DM jacBucketPos = vertcat(std::vector<DM>{horzcat(std::vector<DM>{-lenBA * sin(alpha) - lenAL * sin(alpha + beta) - lenLCoMBucket * sin(alpha + beta + gamma + angMLCoMBucket),
                                                                      -lenAL * sin(alpha + beta) - lenLCoMBucket * sin(alpha + beta + gamma + angMLCoMBucket),
                                                                      -lenLCoMBucket * sin(alpha + beta + gamma + angMLCoMBucket)}),
                                              horzcat(std::vector<DM>{lenBA * cos(alpha) + lenAL * cos(alpha + beta) + lenLCoMBucket * cos(alpha + beta + gamma + angMLCoMBucket),
                                                                      lenAL * cos(alpha + beta) + lenLCoMBucket * cos(alpha + beta + gamma + angMLCoMBucket),
                                                                      lenLCoMBucket * cos(alpha + beta + gamma + angMLCoMBucket)})});
    DM jacBucketRot = horzcat(std::vector<DM>{1, 1, 1});

    DM jacBoomDotPos = vertcat(std::vector<DM>{horzcat(std::vector<DM>{-lenBCoMBoom * cos(alpha + angABCoMBoom) * alphaDot, 0, 0}),
                                               horzcat(std::vector<DM>{-lenBCoMBoom * sin(alpha + angABCoMBoom) * alphaDot, 0, 0})});
    DM jacBoomDotRot = horzcat(std::vector<DM>{0, 0, 0});

    DM jacArmDotPos = vertcat(std::vector<DM>{horzcat(std::vector<DM>{-lenBA * cos(alpha) * alphaDot - lenACoMArm * cos(alpha + beta + angLACoMArm) * (alphaDot + betaDot),
                                                                      -lenACoMArm * cos(alpha + beta + angLACoMArm) * (alphaDot + betaDot),
                                                                      0}),
                                              horzcat(std::vector<DM>{-lenBA * sin(alpha) * alphaDot - lenACoMArm * sin(alpha + beta + angLACoMArm) * (alphaDot + betaDot),
                                                                      -lenACoMArm * sin(alpha + beta + angLACoMArm) * (alphaDot + betaDot),
                                                                      0})});
    DM jacArmDotRot = horzcat(std::vector<DM>{0, 0, 0});

    DM jacBucketDotPos = vertcat(std::vector<DM>{horzcat(std::vector<DM>{-lenBA * cos(alpha) * alphaDot - lenAL * cos(alpha + beta) * (alphaDot + betaDot) - lenLCoMBucket * cos(alpha + beta + gamma + angMLCoMBucket) * (alphaDot + betaDot + gammaDot),
                                                                         -lenAL * cos(alpha + beta) * (alphaDot + betaDot) - lenLCoMBucket * cos(alpha + beta + gamma + angMLCoMBucket) * (alphaDot + betaDot + gammaDot),
                                                                         -lenLCoMBucket * cos(alpha + beta + gamma + angMLCoMBucket) * (alphaDot + betaDot + gammaDot)}),
                                                 horzcat(std::vector<DM>{-lenBA * sin(alpha) * alphaDot - lenAL * sin(alpha + beta) * (alphaDot + betaDot) - lenLCoMBucket * sin(alpha + beta + gamma + angMLCoMBucket) * (alphaDot + betaDot + gammaDot),
                                                                         -lenAL * sin(alpha + beta) * (alphaDot + betaDot) - lenLCoMBucket * sin(alpha + beta + gamma + angMLCoMBucket) * (alphaDot + betaDot + gammaDot),
                                                                         -lenLCoMBucket * sin(alpha + beta + gamma + angMLCoMBucket) * (alphaDot + betaDot + gammaDot)})});
    DM jacBucketDotRot = horzcat(std::vector<DM>{0, 0, 0});

    DM M = mtimes(jacBoomPos.T() * massBoom, jacBoomPos) + mtimes(jacBoomRot.T() * moiBoom, jacBoomRot) +
           mtimes(jacArmPos.T() * massArm, jacArmPos) + mtimes(jacArmRot.T() * moiArm, jacArmRot) +
           mtimes(jacBucketPos.T() * massBucket, jacBucketPos) + mtimes(jacBucketRot.T() * moiBucket, jacBucketRot);
    DM b = mtimes(jacBoomPos.T() * massBoom, mtimes(jacBoomDotPos, qDot)) + mtimes(jacBoomRot.T() * moiBoom, mtimes(jacBoomDotRot, qDot)) +
           mtimes(jacArmPos.T() * massArm, mtimes(jacArmDotPos, qDot)) + mtimes(jacArmRot.T() * moiArm, mtimes(jacArmDotRot, qDot)) +
           mtimes(jacBucketPos.T() * massBucket, mtimes(jacBucketDotPos, qDot)) + mtimes(jacBucketRot.T() * moiBucket, mtimes(jacBucketDotRot, qDot));
    DM g = -mtimes(jacBoomPos.T() * massBoom, gravity) - mtimes(jacArmPos.T() * massArm, gravity) - mtimes(jacBucketPos.T() * massBucket, gravity);

    return mtimes(M, qDDot) + b + g;
}