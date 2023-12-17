#include <vector>
#include "casadi/casadi.hpp"
#include "utils.hpp"

using namespace casadi;

// Converts an angle in radians into a rotation matrix
DM rotMat(float theta)
{
    return vertcat(std::vector<DM>{horzcat(std::vector<DM>{cos(theta), -sin(theta)}),
                                   horzcat(std::vector<DM>{sin(theta), cos(theta)})});
}

// Converts an angle in radians and a position vector into a transformation matrix
DM transMat(float theta, DM pos)
{
    return vertcat(std::vector<DM>{horzcat(std::vector<DM>{rotMat(theta), pos}),
                                   horzcat(std::vector<DM>{0, 0, 1})});
}