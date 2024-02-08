import casadi as csd
import numpy as np
import visualisation as vis
from NLPSolver import NLP, Mode, T, N
from excavatorModel import DutyCycle, inverseKinematics

# n = 10 # Number of time steps
# tSim = 2*T/N

x0 = csd.vertcat(0.0765, -1.8221, -1.2543, 0, 0, 0)
poseDesired = csd.vertcat(1.3, 1, -3)
qDesired = inverseKinematics(poseDesired)

opti = NLP(Mode.NO_LOAD, 0, DutyCycle.S1)
sol = opti.solveNLP(x0, poseDesired)

# qHistory = csd.vertcat(x0[0:3])
# qDotHistory = csd.vertcat(x0[3:6])
# xActual = x0

# for k in range(n):
#     print(k)
#     sol = opti.solveNLP(x0, poseDesired)
#     xActual[0:3] = xActual[0:3] + sol.value(opti.x[3:6, 1])*T/N
#     xActual[3:6] = sol.value(opti.x[3:6, 1])
#     qHistory = csd.horzcat(qHistory, xActual[0:3])
#     qDotHistory = csd.horzcat(qDotHistory, xActual[3:6])
#     x0 = xActual
#     np.savetxt("q.csv", qHistory, delimiter=",")
#     np.savetxt("qDot.csv", qDotHistory, delimiter=',')

# vis.graph(0, tSim, T/N, "Joint Angles", "$\mathsf{q\ (rad)}$", x=sol.value(qHistory))

vis.plotMotorOpPt(sol.value(opti.motorTorque), sol.value(opti.motorVel))
vis.graph(0, T, T/N, "Joint Angles", "$\mathsf{q\ (rad)}$", x=sol.value(opti.x[0:3,:]))
vis.graph(0, T, T/N, "Joint Angular Velocities", "$\mathsf{\dot{q}\ (rad\ s^{-1})}$", x=sol.value(opti.x[3:6,:]))
vis.graph(0, T, T/N, "Joint Angular Accelerations", "$\mathsf{\ddot{q}\ (rad\ s^{-2})}$", u=sol.value(opti.u[:,:]))
vis.graph(0, T, T/N, "Motor Torques", "Motor torque (Nm)", motorTorque=sol.value(opti.motorTorque[:,:]))
vis.graph(0, T, T/N, "Motor Velocities", "Motor velocity ($\mathsf{rad\ s^{-1}}$)", motorVel=sol.value(opti.motorVel[:,:]))
vis.graph(0, T, T/N, "Motor Powers", "Motor power (kW)", motorPower=sol.value(opti.motorPower[:,:])/1000)

for k in range(N + 1):
    if k == 0:
        vis.visualise(sol.value(opti.x[0:3, k]), None, qDesired, k*T/N, k)
    else:
        vis.visualise(sol.value(opti.x[0:3, k]), None, qDesired, k*T/N, k)
    print("Plotting visualisation {k}".format(k=k))

vis.createVideo(0, N, "Excavator", N/T)

print("Done")