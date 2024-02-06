#include <iostream>
#include "casadi/casadi.hpp"
#include "excavatorModel.hpp"
#include "excavatorConstants.hpp"
#include "utils.hpp"

using namespace casadi;

float T = 5; // Time horizon
int N = 50;  // Number of intervals

// System dynamics
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

    DM qInitial = vertcat(DMVector{1, -2.2, -1.8});
    DM poseDesiredValue = vertcat(DMVector{3.37533, -0.897229, -0.4});
    DM qDesired = inverseKinematics(poseDesiredValue);

    opti.set_value(poseDesired, poseDesiredValue);            // Desired pose
    opti.set_value(x0, vertcat(DMVector{qInitial, 0, 0, 0})); // Initial state

    // Weights of loss terms
    MX L = 0;
    float poseAngleFactor = 0.3;
    float angVelFactor = 0.0001;
    float angAccFactor = 0.0001;
    float terminalCostFactor = 1;

    for (int k = 0; k < N; k++)
    {
        // Initial guess
        DM qGuess = (k + 1) / N * (qDesired - qInitial) + qInitial;
        DM qDotGuess = (qDesired - qInitial) / T;

        opti.set_initial(x(all, k), vertcat(DMVector{qGuess, qDotGuess}));

        // Objective function
        MX poseActual = forwardKinematics(static_cast<MX>(x(all, k + 1)));

        L = L + pow(poseActual(0) - poseDesired(0), 2) + pow(poseActual(1) - poseDesired(1), 2) + poseAngleFactor * pow(poseActual(2) - poseDesired(2), 2);
        L = L + angVelFactor * (pow(x(3, k), 2) + pow(x(4, k), 2) + pow(x(5, k), 2));
        L = L + angAccFactor * (pow(u(0, k), 2) + pow(u(1, k), 2) + pow(u(2, k), 2));

        // System dynamics
        opti.subject_to(x(all, k + 1) == integrator(static_cast<MX>(x(all, k)), static_cast<MX>(u(all, k))));

        // Angular constraints
        opti.subject_to(x(Slice(0, 3), k + 1) <= vertcat(DMVector{1.0923, -0.5103, 0.7839}));
        opti.subject_to(x(Slice(0, 3), k + 1) >= vertcat(DMVector{-0.4737, -2.5848, -2.8659}));

        // Motor velocity constraints
        opti.subject_to(motorVel(static_cast<MX>(x(Slice(0, 3), k + 1)),
                                 static_cast<MX>(x(Slice(3, 6), k + 1))) <= vertcat(DMVector{235.62, 235.62, 235.62}));
        opti.subject_to(motorVel(static_cast<MX>(x(Slice(0, 3), k + 1)),
                                 static_cast<MX>(x(Slice(3, 6), k + 1))) >= vertcat(DMVector{-235.62, -235.62, -235.62}));

        // Motor torque constraints
        opti.subject_to(motorTorque(static_cast<MX>(x(Slice(0, 3), k)),
                                    static_cast<MX>(x(Slice(3, 6), k)),
                                    static_cast<MX>(u(all, k)),
                                    vertcat(MXVector{0, 0})) <= vertcat(DMVector{16, 16, 16}));
        opti.subject_to(motorTorque(static_cast<MX>(x(Slice(0, 3), k)),
                                    static_cast<MX>(x(Slice(3, 6), k)),
                                    static_cast<MX>(u(all, k)),
                                    vertcat(MXVector{0, 0})) >= vertcat(DMVector{-16, -16, -16}));
    }

    // Terminal cost
    L = L + terminalCostFactor * (pow(x(3, N), 2) + pow(x(4, N), 2) + pow(x(5, N), 2));

    // Initial state
    opti.subject_to(x(all, 0) == x0);

    opti.minimize(L);
    opti.solver("ipopt");
    OptiSol sol = opti.solve();

    // std::cout << sol.value(x) << std::endl;
    // std::cout << sol.value(u) << std::endl;
}