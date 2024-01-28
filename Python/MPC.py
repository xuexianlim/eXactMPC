import casadi as csd
import excavatorConstants as C
import excavatorModel as mod
import visualisation as vis

T = 0.5 # Time horizon
N = 10 # Number of time steps

def integrator(x, u):
    return csd.vertcat(x[0] + x[3]*T/N + 0.5*u[0]*(T/N)**2,
                       x[1] + x[4]*T/N + 0.5*u[1]*(T/N)**2,
                       x[2] + x[5]*T/N + 0.5*u[2]*(T/N)**2,
                       x[3] + u[0]*T/N,
                       x[4] + u[1]*T/N,
                       x[5] + u[2]*T/N)

opti = csd.Opti()

x = opti.variable(6, N+1)
u = opti.variable(3, N)

__stageLoss__ = opti.variable(N)
#__motorVel__ = opti.variable(3, N)
#__motorTorque__ = opti.variable(3, N)

x0 = opti.parameter(6, 1)
poseDesired = opti.parameter(3, 1)

opti.set_value(poseDesired, csd.vertcat(3.37533, -0.897229, -0.4)) # Desired pose
opti.set_value(x0, csd.vertcat(1, -2.2, -1.8, 0, 0, 0)) # Initial state

# Objective function
L = 0
for k in range(N):
    poseActual = mod.forwardKinematics(x[:, k+1])
    L += (poseActual[0] - poseDesired[0])**2 + (poseActual[1] - poseDesired[1])**2 + (poseActual[2] - poseDesired[2])**2

    opti.subject_to(x[:, k+1] == integrator(x[:, k], u[:, k])) # System dynamics

    opti.subject_to(x[0:3, k+1] <= csd.vertcat(1.0923, -0.5103, 0.7839)) # Angular limits
    opti.subject_to(x[0:3, k+1] >= csd.vertcat(-0.4737, -2.5848, -2.8659))

    #opti.subject_to(x[3:6, k+1] <= csd.vertcat(1, 1, 1)) # Velocity limits
    #opti.subject_to(x[3:6, k+1] >= csd.vertcat(-1, -1, -1))

    opti.subject_to(u[:, k] <= csd.vertcat(20, 20, 20)) # Acceleration limits
    opti.subject_to(u[:, k] >= csd.vertcat(-20, -20, -20))

    #opti.subject_to(motorVel(x[0:3, k], x[3:6, k]) <= csd.vertcat(235.62, 235.62, 235.62)) # Motor speed
    #opti.subject_to(motorVel(x[0:3, k], x[3:6, k]) >= csd.vertcat(-235.62, -235.62, -235.62)) # 2250 RPM

    #opti.subject_to(motorTorque(x[0:3, k], x[3:6, k], u[:, k]) <= csd.vertcat(16, 16, 16)) # Motor torque
    #opti.subject_to(motorTorque(x[0:3, k], x[3:6, k], u[:, k]) >= csd.vertcat(-16, -16, -16))

    # For printing and viewing variables of interest
    opti.subject_to(__stageLoss__[k] ==  (poseActual[0] - poseDesired[0])**2 + (poseActual[1] - poseDesired[1])**2 + (poseActual[2] - poseDesired[2])**2)
    #opti.subject_to(__motorVel__[:, k] == mod.motorVel(x[0:3, k], x[3:6, k]))
    #opti.subject_to(__motorTorque__[:, k] == mod.motorTorque(x[0:3, k], x[3:6, k], u[:, k]))

opti.subject_to(x[:, 0] == x0) # Initial state

opts={} 
# opts["verbose_init"] = False
# opts["verbose"] = False
# opts["print_time"] = False
# opts["ipopt.print_level"] = 0

opti.minimize(L)
opti.solver('ipopt', opts)

sol = opti.solve()

print("==============================   x   ==============================")
print(sol.value(x))
print("==============================   u   ==============================")
print(sol.value(u))

vis.graph(0, T, N, "Test", x=sol.value(x), u=sol.value(u))
vis.graph(0, T, N, "Loss", stageLoss=sol.value(__stageLoss__))

for k in range(N + 1):
    if k == 0:
        vis.visualise(sol.value(x[0:3, k]), None, csd.vertcat(0, -0.7, 0.3), k*T/N, k)
    else:
        vis.visualise(sol.value(x[0:3, k]), sol.value(x[0:3, k - 1]), csd.vertcat(0, -0.7, 0.3), k*T/N, k)

print("Done")