#ifndef EXCAVATOR_CONSTANTS_H
#define EXCAVATOR_CONSTANTS_H

#include "casadi/casadi.hpp"

using namespace casadi;

// Boom
const DM iBC = vertcat(DMVector{0.135, -0.264});
const DM bBA = vertcat(DMVector{2.050, 0});
const DM bBD = vertcat(DMVector{1.025, 0.278});
const DM bBE = vertcat(DMVector{0.977, 0.737});
const DM bDA = bBA - bBD;
const DM bAE = bBE - bBA;
const DM bCoMBoom = vertcat(DMVector{1.025, 0.384});

const DM lenBC = norm_2(iBC);
const DM lenBA = norm_2(bBA);
const DM lenBD = norm_2(bBD);
const DM lenBE = norm_2(bBE);
const DM lenDA = norm_2(bDA);
const DM lenAE = norm_2(bAE);
const DM lenBCoMBoom = norm_2(bCoMBoom);

const DM angABD = acos((pow(lenBA, 2) + pow(lenBD, 2) - pow(lenDA, 2)) /
                       (2 * lenBA * lenBD));
const DM angABCoMBoom = atan2(bCoMBoom(1), bCoMBoom(0));

const float massBoom = 227.343;
const float moiBoom = 67.768;

// Arm
const DM aAF = vertcat(DMVector{-0.251, 0.158});
const DM aAG = vertcat(DMVector{-0.134, 0.320});
const DM aAJ = vertcat(DMVector{0.880, 0});
const DM aAL = vertcat(DMVector{1.050, 0});
const DM aFL = aAL - aAF;
const DM aJG = aAG - aAJ;
const DM aJL = aAL - aAJ;
const DM aCoMArm = vertcat(DMVector{0.225, 0.227});

const float lenHJ = 0.240;
const float lenHK = 0.240;
const DM lenAF = norm_2(aAF);
const DM lenAG = norm_2(aAG);
const DM lenAJ = norm_2(aAJ);
const DM lenAL = norm_2(aAL);
const DM lenFL = norm_2(aFL);
const DM lenJG = norm_2(aJG);
const DM lenJL = norm_2(aJL);
const DM lenACoMArm = norm_2(aCoMArm);

const DM angFAL = acos((pow(lenAL, 2) + pow(lenAF, 2) - pow(lenFL, 2)) /
                       (2 * lenAF * lenAL));
const DM angLACoMArm = atan2(aCoMArm(1), aCoMArm(0));

const float massArm = 130.123;
const float moiArm = 30.258;

// Bucket
const DM lLK = vertcat(DMVector{-0.014, 0.164});
const DM lLM = vertcat(DMVector{0.567, 0});
const DM lKM = lLM - lLK;
const DM lCoMBucket = vertcat(DMVector{0.289, 0.166});

const DM lenLK = norm_2(lLK);
const DM lenLM = norm_2(lLM);
const DM lenKM = norm_2(lKM);
const DM lenLCoMBucket = norm_2(lCoMBucket);

const DM angKLM = acos((pow(lenLM, 2) + pow(lenLK, 2) - pow(lenKM, 2)) /
                       (2 * lenLK * lenLM));
const DM angMLCoMBucket = atan2(lCoMBucket(1), lCoMBucket(0));

const float massBucket = 53.000;
const float moiBucket = 3.021;

// Environment
const float yGround = -0.95737;
const DM gravity = vertcat(DMVector{0, -9.81});

#endif