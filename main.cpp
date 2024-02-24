#include <iostream>
#include "casadi/casadi.hpp"
#include "excavatorModel.hpp"
#include "excavatorConstants.hpp"
#include "utils.hpp"

using namespace casadi;

float T = 2;      // Time horizon
int N = 20;       // Number of intervals
float Ts = T / N; // MPC sampling time

Mode mode = LIFT;
int extF = 1000;
DutyCycle dutyCycle = S2_30;

// System dynamics
template <typename V>
V integrator(V x, V u, float t)
{
    return vertcat(std::vector<V>{x(0) + x(3) * t + 0.5 * u(0) * pow(t, 2),
                                  x(1) + x(4) * t + 0.5 * u(1) * pow(t, 2),
                                  x(2) + x(5) * t + 0.5 * u(2) * pow(t, 2),
                                  x(3) + u(0) * t,
                                  x(4) + u(1) * t,
                                  x(5) + u(2) * t});
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

    // Pseudo "states" for convenience
    MX __motorVel__ = opti.variable(3, N);
    MX __motorTorque__ = opti.variable(3, N);
    MX __motorPower__ = opti.variable(3, N);

    DM qInitial = vertcat(DMVector{1, -2.2, -1.8});
    DM poseDesiredValue = vertcat(DMVector{3.37533, -0.897229, -0.4});
    DM qDesired = inverseKinematics(poseDesiredValue);
    DM extForce = extF;

    opti.set_value(poseDesired, poseDesiredValue);            // Desired pose
    opti.set_value(x0, vertcat(DMVector{qInitial, 0, 0, 0})); // Initial state

    // Weights of loss terms
    float xWeight = 1;
    float yWeight = 1;
    float thetaWeight = 0.3;
    float terminalCostWeight = 1;
    float regularisationWeight = 0.01;

    // Objective function
    MX L = 0;

    for (int k = 0; k < N; k++)
    {
        // Initial guess
        DM qGuess = (k + 1) / N * (qDesired - qInitial) + qInitial;
        DM qDotGuess = (qDesired - qInitial) / T;

        opti.set_initial(x(all, k), vertcat(DMVector{qGuess, qDotGuess}));

        // Objective function
        MX poseActual = forwardKinematics(static_cast<MX>(x(all, k + 1)));

        L = L + xWeight * pow(poseActual(0) - poseDesired(0), 2) + yWeight * pow(poseActual(1) - poseDesired(1), 2) + thetaWeight * pow(poseActual(2) - poseDesired(2), 2);
        L = L + regularisationWeight * (pow(u(0, k), 2) + pow(u(1, k), 2) + pow(u(2, k), 2));

        // System dynamics
        opti.subject_to(x(all, k + 1) == integrator(static_cast<MX>(x(all, k)), static_cast<MX>(u(all, k)), Ts));
        opti.subject_to(__motorVel__(all, k) == motorVel(static_cast<MX>(x(Slice(0, 3), k)), static_cast<MX>(x(Slice(3, 6), k))));
        if (mode == LIFT)
        {
            opti.subject_to(__motorTorque__(all, k) == motorTorque(static_cast<MX>(x(Slice(0, 3), k)), static_cast<MX>(x(Slice(3, 6), k)), static_cast<MX>(u(all, k)), vertcat(MXVector{0, -extForce})));
        }
        else
        {
            opti.subject_to(__motorTorque__(all, k) == motorTorque(static_cast<MX>(x(Slice(0, 3), k)), static_cast<MX>(x(Slice(3, 6), k)), static_cast<MX>(u(all, k)), vertcat(MXVector{0, 0})));
        }
        opti.subject_to(__motorPower__(all, k) == __motorTorque__(all, k) * __motorVel__(all, k));

        // Angular constraints
        opti.subject_to(x(Slice(0, 3), k + 1) <= vertcat(DMVector{1.0923, -0.5103, 0.7839}));
        opti.subject_to(x(Slice(0, 3), k + 1) >= vertcat(DMVector{-0.4737, -2.5848, -2.8659}));

        // Motor velocity constraints
        opti.subject_to(__motorVel__(all, k) <= vertcat(DMVector{471.24, 471.24, 471.24}));
        opti.subject_to(__motorVel__(all, k) >= vertcat(DMVector{-471.24, -471.24, -471.24}));

        // Motor torque constraints
        opti.subject_to(__motorTorque__(all, k) <= motorTorqueLimit(static_cast<MX>(__motorVel__(all, k)), dutyCycle));
        opti.subject_to(__motorTorque__(all, k) >= -motorTorqueLimit(static_cast<MX>(__motorVel__(all, k)), dutyCycle));
        opti.subject_to(__motorTorque__(all, k) <= motorTorqueLimit(static_cast<MX>(-__motorVel__(all, k)), dutyCycle));
        opti.subject_to(__motorTorque__(all, k) >= -motorTorqueLimit(static_cast<MX>(-__motorVel__(all, k)), dutyCycle));

        // Motor power constraints
        opti.subject_to(__motorPower__(0, k) + __motorPower__(1, k) + __motorPower__(2, k) <= 8250);
        opti.subject_to(__motorPower__(0, k) + __motorPower__(1, k) - __motorPower__(2, k) <= 8250);
        opti.subject_to(__motorPower__(0, k) - __motorPower__(1, k) + __motorPower__(2, k) <= 8250);
        opti.subject_to(__motorPower__(0, k) - __motorPower__(1, k) - __motorPower__(2, k) <= 8250);
        opti.subject_to(-__motorPower__(0, k) + __motorPower__(1, k) + __motorPower__(2, k) <= 8250);
        opti.subject_to(-__motorPower__(0, k) + __motorPower__(1, k) - __motorPower__(2, k) <= 8250);
        opti.subject_to(-__motorPower__(0, k) - __motorPower__(1, k) + __motorPower__(2, k) <= 8250);
        opti.subject_to(-__motorPower__(0, k) - __motorPower__(1, k) - __motorPower__(2, k) <= 8250);
    }

    // Terminal cost
    L = L + terminalCostWeight * (pow(x(3, N), 2) + pow(x(4, N), 2) + pow(x(5, N), 2));

    // Initial state
    opti.subject_to(x(all, 0) == x0);

    opti.minimize(L);
    opti.solver("ipopt");
    OptiSol sol = opti.solve();

    float TMotor = 0.001;     // Motor velocity command time interval
    int NMotor = Ts / TMotor; // Number of motor velocity commands per MPC time step

    for (int i = 0; i < NMotor; i++)
    {
        DM xInterpolate = integrator(sol.value(x(all, 0)), sol.value(u(all, 0)), (i + 1) * TMotor);
        DM motorSpd = motorVel(static_cast<DM>(xInterpolate(Slice(0, 3))), static_cast<DM>((xInterpolate(Slice(0, 3)))));

        std::cout << motorSpd << std::endl; // Motor velocity commands every TMotor seconds
    }
}