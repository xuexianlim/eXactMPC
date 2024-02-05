#include "utils.hpp"

using namespace casadi;

// Converts an angle in radians into a rotation matrix
DM rotMat(float theta)
{
    return vertcat(DMVector{horzcat(DMVector{cos(theta), -sin(theta)}),
                            horzcat(DMVector{sin(theta), cos(theta)})});
}

// Converts an angle in radians and a position vector into a transformation matrix
DM transMat(float theta, DM pos)
{
    return vertcat(DMVector{horzcat(DMVector{rotMat(theta), pos}),
                            horzcat(DMVector{0, 0, 1})});
}