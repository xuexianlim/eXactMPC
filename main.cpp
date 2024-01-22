#include <iostream>
#include <vector>
#include "casadi/casadi.hpp"
#include "excavatorModel.hpp"
#include "excavatorConstants.hpp"
#include "utils.hpp"

using namespace casadi;

float T = 5; // Time horizon
int N = 20;  // Number of intervals

template <typename V>
V integrator(V x, V u)
{
    return vertcat(std::vector<V>{x(0) + x(3) * T / N + 0.5 * u(0) * pow(T / N, 2),
                                  x(1) + x(4) * T / N + 0.5 * u(1) * pow(T / N, 2),
                                  x(2) + x(5) * T / N + 0.5 * u(2) * pow(T / N, 2),
                                  x(3) + u(0) * T / N,
                                  x(4) + u(1) * T / N,
                                  x(5) + u(2) * T / N});
}

int main()
{
    Opti opti = Opti(); // Create optimisation problem

    Slice all;

    // Optimisation variables
    MX x = opti.variable(6, N + 1);
    MX u = opti.variable(3, N);
    MX x0 = opti.parameter(6, 1);
    MX poseDesired = opti.parameter(3, 1);

    opti.set_value(poseDesired, vertcat(DMVector{0, 3.667, pi / 2}));
    opti.set_value(x0, vertcat(DMVector{0, -1, 0, 0, 0, 0}));

    MX L = 0;
    for (int k = 0; k < N; k++)
    {
        // Objective function
        MX poseActual = forwardKinematics(static_cast<MX>(x(all, k + 1)));

        L = L + pow(poseActual(0) - poseDesired(0), 2) + pow(poseActual(1) - poseDesired(1), 2) + pow(poseActual(2) - poseDesired(2), 2);

        // System dynamics
        opti.subject_to(x(all, k + 1) == integrator(static_cast<MX>(x(all, k)), static_cast<MX>(u(all, k))));

        // Angular constraints
        opti.subject_to(x(Slice(0, 3), k + 1) <= vertcat(DMVector{1.0923, -0.5103, 0.7839}));
        opti.subject_to(x(Slice(0, 3), k + 1) >= vertcat(DMVector{-0.4737, -2.5848, -2.8659}));

        // Angular velocity constraints
        opti.subject_to(x(Slice(3, 6), k + 1) <= vertcat(DMVector{1, 1, 1}));
        opti.subject_to(x(Slice(3, 6), k + 1) >= vertcat(DMVector{-1, -1, -1}));

        // Angular acceleration constraints
        opti.subject_to(u(all, k) <= vertcat(DMVector{1, 1, 1}));
        opti.subject_to(u(all, k) >= vertcat(DMVector{-1, -1, -1}));

        // Torque constraints
        opti.subject_to(inverseDynamics(static_cast<MX>(x(Slice(0, 3), k)),
                                        static_cast<MX>(x(Slice(3, 6), k)),
                                        static_cast<MX>(u(all, k))) <= vertcat(DMVector{30000, 30000, 30000}));
    }
    // Initial state
    opti.subject_to(x(all, 0) == x0);

    opti.minimize(L);
    opti.solver("ipopt");
    OptiSol sol = opti.solve();
    std::cout << sol.value(x) << std::endl;

    return 0;
}
