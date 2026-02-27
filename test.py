import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider
import matplotlib
matplotlib.use("TkAgg")

class SimpleMotor:
    def __init__(self, dt=0.001, duration=2.0):
        self.I_arm = 0.01        
        self.I_motor = 0.001     
        self.gear_ratio = 10.0 
        self.damping = 0.1      
        
        self.kp = 50.0        
        self.kd = 5.0           
        self.tau_max = 5.0      
        
        self.dt = dt
        self.steps = int(duration / dt)
        self.time = np.linspace(0, duration, self.steps)
        
        self.I_eff = self.I_arm + (self.I_motor * (self.gear_ratio ** 2))

    def simulate(self, theta_target):

        theta = 0.0      
        theta_dot = 0.0  
        
        theta_hist = np.zeros(self.steps)
        theta_dot_hist = np.zeros(self.steps)
        theta_ddot_hist = np.zeros(self.steps)
        torque_hist = np.zeros(self.steps)

        for i in range(self.steps):

            error = theta_target - theta
            
            motor_torque = (self.kp * error) - (self.kd * theta_dot)

            motor_torque = np.clip(motor_torque, -self.tau_max, self.tau_max)

            joint_torque = motor_torque * self.gear_ratio
            theta_ddot = (joint_torque - (self.damping * theta_dot)) / self.I_eff

            theta_dot += theta_ddot * self.dt
            theta += theta_dot * self.dt

            theta_hist[i] = theta
            theta_dot_hist[i] = theta_dot
            theta_ddot_hist[i] = theta_ddot
            torque_hist[i] = motor_torque

        return theta_hist, theta_dot_hist, theta_ddot_hist, torque_hist

motor = SimpleMotor()
theta_target = np.pi / 2 
theta_hist, theta_dot_hist, theta_ddot_hist, torque_hist = motor.simulate(theta_target)

fig = plt.figure(figsize=(12, 9))
gs = fig.add_gridspec(4, 2)

ax_arm = fig.add_subplot(gs[:, 0])
ax_pos = fig.add_subplot(gs[0, 1])
ax_vel = fig.add_subplot(gs[1, 1])
ax_acc = fig.add_subplot(gs[2, 1])
ax_tau = fig.add_subplot(gs[3, 1])

ax_arm.set_aspect("equal")
ax_arm.set_xlim(-0.4, 0.4); ax_arm.set_ylim(-0.4, 0.4)
arm_line, = ax_arm.plot([], [], lw=5, color='blue', solid_capstyle='round')
ax_arm.add_patch(plt.Circle((0, 0), 0.05, color='gray'))


pos_line, = ax_pos.plot([], [], 'r', label="θ [rad]")
target_line = ax_pos.axhline(theta_target, ls='--', color='black', alpha=0.5)
ax_pos.set_xlim(0, 2); ax_pos.set_ylim(-0.2, 4.0)
ax_pos.set_title("Vị trí (Position)")

vel_line, = ax_vel.plot([], [], 'g')
ax_vel.set_xlim(0, 2); ax_vel.set_ylim(-20, 20)
ax_vel.set_title("Vận tốc (Velocity)")

acc_line, = ax_acc.plot([], [], 'm')
ax_acc.set_xlim(0, 2); ax_acc.set_ylim(-100, 100)
ax_acc.set_title("Gia tốc (Acceleration)")

tau_line, = ax_tau.plot([], [], 'orange')
ax_tau.set_xlim(0, 2); ax_tau.set_ylim(-6, 6)
ax_tau.set_title("Mô-men xoắn (Motor Torque)")

def update(frame):

    angle = theta_hist[frame]
    arm_line.set_data([0, 0.3 * np.cos(angle)], [0, 0.3 * np.sin(angle)])
    

    t = motor.time[:frame]
    pos_line.set_data(t, theta_hist[:frame])
    vel_line.set_data(t, theta_dot_hist[:frame])
    acc_line.set_data(t, theta_ddot_hist[:frame])
    tau_line.set_data(t, torque_hist[:frame])
    return arm_line, pos_line, vel_line, acc_line, tau_line

ani = FuncAnimation(fig, update, frames=range(0, motor.steps, 10), interval=10, blit=True)

ax_slider = plt.axes([0.15, 0.02, 0.7, 0.03])
slider = Slider(ax_slider, 'Target (rad)', 0, np.pi, valinit=theta_target)

def change_target(val):
    global theta_hist, theta_dot_hist, theta_ddot_hist, torque_hist
    target_line.set_ydata([val, val])
    theta_hist, theta_dot_hist, theta_ddot_hist, torque_hist = motor.simulate(val)
    ani.frame_seq = ani.new_frame_seq()

slider.on_changed(change_target)
plt.subplots_adjust(bottom=0.15)
plt.show()
