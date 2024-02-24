import casadi as csd
from enum import Enum
import excavatorConstants as C
import excavatorModel as mod

T = 2 # Prediction horizon
N = 20 # Number of time steps
Ts = T/N # Sampling time

# Changes how the external force behaves
class Mode(Enum):
    NO_LOAD = 0 # No external forces
    LIFT = 1 # External load in the bucket exerts a downward force due to gravity always
    DIG = 2 # Resistive force acts against the direction of the velocity of the bucket

def integrator(x, u, t):
    return csd.vertcat(x[0] + x[3]*t + 0.5*u[0]*(t)**2,
                       x[1] + x[4]*t + 0.5*u[1]*(t)**2,
                       x[2] + x[5]*t + 0.5*u[2]*(t)**2,
                       x[3] + u[0]*t,
                       x[4] + u[1]*t,
                       x[5] + u[2]*t)

class NLP():
    def __init__(self, mode, extF, dutyCycle):
        self.opti = csd.Opti()

        self.x = self.opti.variable(6, N+1) # State
        self.u = self.opti.variable(3, N) # Input

        self.x0 = self.opti.parameter(6, 1) # Initial state
        self.poseDesired = self.opti.parameter(3, 1) # Desired pose
        self.extForce = self.opti.parameter() # External force magnitude

        # Other "states" for convenience and to avoid repetition of the same calculations by the solver
        self.poseActual = self.opti.variable(3, N)
        self.extForce = self.opti.variable(2, N)
        self.motorTorque = self.opti.variable(3, N)
        self.motorVel = self.opti.variable(3, N)
        self.motorPower = self.opti.variable(3, N)

        # Objective function weights
        xWeight = 1
        yWeight = 1
        thetaWeight = C.lenLM**2
        terminalWeight = 1
        regularisationWeight = 0.01

        # Objective function
        L = 0

        for k in range(N):
            # Objective function
            L += xWeight*(self.poseActual[0, k] - self.poseDesired[0])**2 + yWeight*(self.poseActual[1, k] - self.poseDesired[1])**2 + thetaWeight*(self.poseActual[2, k] - self.poseDesired[2])**2 # Cost due to difference from desired pose
            L += regularisationWeight*(self.u[0, k]**2 + self.u[1, k]**2 + self.u[2, k]**2) # Regularisation cost

            # Constraints
            self.opti.subject_to(self.x[:, k+1] == integrator(self.x[:, k], self.u[:, k], Ts)) # System dynamics
            self.opti.subject_to(self.poseActual[:, k] == mod.forwardKinematics(self.x[:, k+1]))
            self.opti.subject_to(self.motorTorque[:, k] == mod.motorTorque(self.x[0:3, k], self.x[3:6, k], self.u[:, k], self.extForce))
            self.opti.subject_to(self.motorVel[:, k] == mod.motorVel(self.x[0:3, k], self.x[3:6, k]))
            self.opti.subject_to(self.motorPower[:, k] == self.motorTorque[:, k]*self.motorVel[:, k])
            if mode == Mode.LIFT:
                self.opti.subject_to(self.extForce[:, k] == csd.vertcat(0, -extF))
            #elif mode == Mode.DIG:
                #Not implemented
            else:
                self.opti.subject_to(self.extForce[:, k] == csd.vertcat(0, 0))

            self.opti.subject_to(self.x[0:3, k+1] <= csd.vertcat(1.0923, -0.5103, 0.7839)) # Angular limits
            self.opti.subject_to(self.x[0:3, k+1] >= csd.vertcat(-0.4737, -2.5848, -2.8659))

            self.opti.subject_to(self.motorVel[:, k] <= csd.vertcat(471.24, 471.24, 471.24)) # Motor speed
            self.opti.subject_to(self.motorVel[:, k] >= csd.vertcat(-471.24, -471.24, -471.24)) # 4500 RPM

            self.opti.subject_to(self.motorTorque[:, k] <= mod.motorTorqueLimit(self.motorVel[:, k], dutyCycle)) # Motor torque limits
            self.opti.subject_to(self.motorTorque[:, k] >= -mod.motorTorqueLimit(self.motorVel[:, k], dutyCycle))
            self.opti.subject_to(self.motorTorque[:, k] <= mod.motorTorqueLimit(-self.motorVel[:, k], dutyCycle))
            self.opti.subject_to(self.motorTorque[:, k] >= -mod.motorTorqueLimit(-self.motorVel[:, k], dutyCycle))

            motorPower = 8250
            self.opti.subject_to(self.motorPower[0, k] + self.motorPower[1, k] + self.motorPower[2, k] <= motorPower) # Motor power limits
            self.opti.subject_to(self.motorPower[0, k] + self.motorPower[1, k] - self.motorPower[2, k] <= motorPower)
            self.opti.subject_to(self.motorPower[0, k] - self.motorPower[1, k] + self.motorPower[2, k] <= motorPower)
            self.opti.subject_to(self.motorPower[0, k] - self.motorPower[1, k] - self.motorPower[2, k] <= motorPower)
            self.opti.subject_to(-self.motorPower[0, k] + self.motorPower[1, k] + self.motorPower[2, k] <= motorPower)
            self.opti.subject_to(-self.motorPower[0, k] + self.motorPower[1, k] - self.motorPower[2, k] <= motorPower)
            self.opti.subject_to(-self.motorPower[0, k] - self.motorPower[1, k] + self.motorPower[2, k] <= motorPower)
            self.opti.subject_to(-self.motorPower[0, k] - self.motorPower[1, k] - self.motorPower[2, k] <= motorPower)

        L += terminalWeight*(self.x[3, N]**2 + self.x[4, N]**2 + self.x[5, N]**2)

        self.opti.minimize(L)

    def solveNLP(self, x0, poseDesired):
        # Initial state and desired pose
        self.opti.set_value(self.x0, x0)
        self.opti.set_value(self.poseDesired, poseDesired)

        self.opti.subject_to(self.x[:, 0] == self.x0)

        # Initial guess for x
        for k in range(N):
            qDesired = mod.inverseKinematics(poseDesired)
            qGuess = (k+1)/N*(qDesired - x0[0:3]) + x0[0:3]
            qDotGuess = (qDesired - x0[0:3])/T

            self.opti.set_initial(self.x[:, k], csd.vertcat(qGuess, qDotGuess))

        # Options
        opts={}

        # Comment out this block to see Ipopt solver printouts
        opts["verbose_init"] = False
        opts["verbose"] = False
        opts["print_time"] = False
        opts["ipopt.print_level"] = 0

        self.opti.solver('ipopt', opts)
        sol = self.opti.solve()

        return sol