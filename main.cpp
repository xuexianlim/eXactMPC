#include <iostream>
#include <vector>
#include "casadi/casadi.hpp"
#include "excavatorModel.hpp"
#include "excavatorConstants.hpp"
#include "utils.hpp"

using namespace casadi;

int main() {

    DM q = vertcat(std::vector<DM>{1.0923,-0.5103,-2.8659});
    DM qDot = vertcat(std::vector<DM>{1,1,1});
    DM qDDot = vertcat(std::vector<DM>{1,1,1});
    std::cout << inverseDynamics(q, qDot, qDDot) << std::endl;
    
    return 0;
}