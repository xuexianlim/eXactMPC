#ifndef UTILS_H
#define UTILS_H

#include "casadi/casadi.hpp"

using namespace casadi;

DM rotMat(float theta); 
DM transMat(float theta, DM pos);

#endif