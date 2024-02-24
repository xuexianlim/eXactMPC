import casadi as csd
import visualisation as vis
from NLPSolver import NLP, Mode, T, N
from excavatorModel import DutyCycle, inverseKinematics

x0 = csd.vertcat(1, -2.2, -1.8, 0, 0, 0)
poseDesired = csd.vertcat(3.2, -0.95737, -0.8)
qDesired = inverseKinematics(poseDesired)

mode = Mode.LIFT
dutyCycle = DutyCycle.S2_30

opti = NLP(mode, 500, dutyCycle)
sol = opti.solveNLP(x0, poseDesired)

vis.plotMotorOpPt(sol.value(opti.motorTorque), sol.value(opti.motorVel), dutyCycle)
vis.graph(0, T, T/N, "Joint Angles", "$\mathsf{q\ (rad)}$", x=sol.value(opti.x[0:3,:]))
vis.graph(0, T, T/N, "Joint Angular Velocities", "$\mathsf{\dot{q}\ (rad\ s^{-1})}$", x=sol.value(opti.x[3:6,:]))
vis.graph(0, T, T/N, "Joint Angular Accelerations", "$\mathsf{\ddot{q}\ (rad\ s^{-2})}$", u=sol.value(opti.u[:,:]))
vis.graph(0, T, T/N, "Motor Torques", "Motor torque (Nm)", motorTorque=sol.value(opti.motorTorque[:,:]))
vis.graph(0, T, T/N, "Motor Velocities", "Motor velocity ($\mathsf{rad\ s^{-1}}$)", motorVel=sol.value(opti.motorVel[:,:]))
vis.graph(0, T, T/N, "Motor Powers", "Motor power (kW)", motorPower=sol.value(opti.motorPower[:,:])/1000)

for k in range(N + 1):
    if k == 0:
        if mode == Mode.LIFT:
            vis.visualise(sol.value(opti.x[0:3, k]), None, qDesired, k*T/N, k, sol.value(opti.extForce[:, 0]))
        else:
            vis.visualise(sol.value(opti.x[0:3, k]), None, qDesired, k*T/N, k, None)
    else:
        if mode == Mode.LIFT:
            vis.visualise(sol.value(opti.x[0:3, k]), None, qDesired, k*T/N, k, sol.value(opti.extForce[:, k-1]))
        else:
            vis.visualise(sol.value(opti.x[0:3, k]), None, qDesired, k*T/N, k, None)
    print("Plotting visualisation {k}".format(k=k))

vis.createVideo(0, N, "Excavator", N/T)

print("Done")