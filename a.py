"""
Chapter 2: Simple Oscillator Simulation
========================================

A simplified mass-spring-damper system demonstrating:
- Real-time animation with spring visualization
- Position and velocity plots over time
- Interactive parameter adjustment (mass, damping)

Physical System: m·ẍ + d·ẋ + k·x = 0
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider, Button
from matplotlib.patches import Rectangle


# =====================================================
# Simple Oscillator Model
# =====================================================
class SimpleOscillator:
    """
    Mass-spring-damper system (1D oscillation)

    Technical Variables:
        - position [m]: displacement from equilibrium
        - velocity [m/s]: rate of change of position
        - acceleration [m/s²]: rate of change of velocity
    """

    def __init__(self, mass, spring_constant, damping, initial_position):
        """Initialize oscillator with parameters."""

        # Validate parameters
        if mass <= 0:
            raise ValueError("Mass must be positive.")
        if spring_constant <= 0:
            raise ValueError("Spring constant must be positive.")
        if damping < 0:
            raise ValueError("Damping cannot be negative.")

        # System parameters
        self.mass = mass
        self.k = spring_constant
        self.d = damping

        # State variables
        self.position = initial_position
        self.velocity = 0.0
        self.time = 0.0

        # Store initial conditions
        self.initial_position = initial_position
        self.initial_velocity = 0.0

    def compute_system_info(self):
        """Calculate system parameters for display."""

        omega_0 = np.sqrt(self.k / self.mass)
        damping_ratio = self.d / (2 * np.sqrt(self.k * self.mass))

        if damping_ratio == 0:
            damping_type = "Undamped"
        elif damping_ratio < 1:
            damping_type = "Underdamped"
        elif np.isclose(damping_ratio, 1.0):
            damping_type = "Critically damped"
        else:
            damping_type = "Overdamped"

        return {
            'omega_0': omega_0,
            'damping_ratio': damping_ratio,
            'damping_type': damping_type
        }

    def compute_acceleration(self):
        """Calculate acceleration using Newton's law."""

        acceleration = -(self.d / self.mass) * self.velocity \
                       - (self.k / self.mass) * self.position

        return acceleration

    def step(self, dt):
        """Advance simulation using Euler integration."""

        accel = self.compute_acceleration()

        # Euler integration
        self.velocity = self.velocity + accel * dt
        self.position = self.position + self.velocity * dt
        self.time = self.time + dt

        return self.position, self.velocity, accel

    def reset(self):
        """Reset oscillator to initial conditions."""

        self.position = self.initial_position
        self.velocity = self.initial_velocity
        self.time = 0.0


# =====================================================
# Visualization
# =====================================================
class OscillatorVisualizer:
    """Interactive visualization for the oscillator."""

    def __init__(self, oscillator, dt=0.002, total_time=15.0):
        self.osc = oscillator
        self.dt = dt
        self.total_time = total_time
        self.steps = int(total_time / dt)
        self.time_array = np.linspace(0, total_time, self.steps)

        self.position_history = np.zeros(self.steps)
        self.velocity_history = np.zeros(self.steps)
        self.frame = 0

        self.setup_figure()
        self.create_controls()
        self.update_info_display()

    def setup_figure(self):
        self.fig = plt.figure(figsize=(12, 8))
        self.fig.suptitle('Simple Oscillator Simulation - Chapter 2',
                          fontsize=14, fontweight='bold')

        gs = self.fig.add_gridspec(3, 1, height_ratios=[1.5, 1, 1],
                                   left=0.1, right=0.9, top=0.92, bottom=0.25,
                                   hspace=0.4)

        self.ax_anim = self.fig.add_subplot(gs[0])
        self.setup_animation()

        self.ax_pos = self.fig.add_subplot(gs[1])
        self.setup_position_plot()

        self.ax_vel = self.fig.add_subplot(gs[2])
        self.setup_velocity_plot()

    def setup_animation(self):
        self.ax_anim.set_xlim(-0.15, 0.15)
        self.ax_anim.set_ylim(-0.05, 0.15)
        self.ax_anim.set_title("Mass-Spring-Damper Animation",
                               fontsize=12, fontweight='bold')
        self.ax_anim.set_aspect('equal')
        self.ax_anim.axis('off')

        self.ax_anim.plot([-0.12, -0.12], [0.08, 0.12], 'k-', linewidth=10)

        self.spring_line, = self.ax_anim.plot([], [], 'b-', linewidth=2.5)

        self.mass_rect = Rectangle((0, 0), 0.035, 0.035,
                                   facecolor='steelblue',
                                   edgecolor='darkblue',
                                   linewidth=2.5)
        self.ax_anim.add_patch(self.mass_rect)

        self.damper_line, = self.ax_anim.plot([], [], 'r-', linewidth=2)

        self.ax_anim.axvline(x=0, color='green', linestyle='--',
                             linewidth=1.5, alpha=0.4)

        self.pos_text = self.ax_anim.text(0.0, 0.02, '', fontsize=11)

        self.info_text = self.ax_anim.text(-0.13, -0.03, '', fontsize=9,
                                           family='monospace',
                                           verticalalignment='top')

    def draw_spring(self, x_start, x_end, y=0.1, n_coils=10):
        length = x_end - x_start
        coil_amplitude = 0.012

        n_points = n_coils * 4 + 2
        spring_x = np.zeros(n_points)
        spring_y = np.zeros(n_points)

        spring_x[0] = x_start
        spring_y[0] = y

        for i in range(1, n_points - 1):
            progress = i / (n_points - 1)
            spring_x[i] = x_start + progress * length

            phase = (i - 1) % 4
            if phase == 1:
                spring_y[i] = y + coil_amplitude
            elif phase == 3:
                spring_y[i] = y - coil_amplitude
            else:
                spring_y[i] = y

        spring_x[-1] = x_end
        spring_y[-1] = y

        return spring_x, spring_y

    def setup_position_plot(self):
        self.ax_pos.set_xlim(0, self.total_time)
        self.ax_pos.set_ylim(-0.1, 0.1)
        self.ax_pos.set_title('Position vs Time')
        self.ax_pos.grid(True)

        self.pos_line, = self.ax_pos.plot([], [], 'b-')

    def setup_velocity_plot(self):
        self.ax_vel.set_xlim(0, self.total_time)
        self.ax_vel.set_ylim(-1.5, 1.5)
        self.ax_vel.set_title('Velocity vs Time')
        self.ax_vel.grid(True)

        self.vel_line, = self.ax_vel.plot([], [], 'g-')

    def create_controls(self):
        ax_mass = plt.axes([0.2, 0.14, 0.6, 0.02])
        self.slider_mass = Slider(ax_mass, 'Mass [kg]',
                                  0.1, 2.0,
                                  valinit=self.osc.mass,
                                  valstep=0.05)

        ax_damping = plt.axes([0.2, 0.10, 0.6, 0.02])
        self.slider_damping = Slider(ax_damping, 'Damping [Ns/m]',
                                     0.0, 5.0,
                                     valinit=self.osc.d,
                                     valstep=0.1)

        ax_x0 = plt.axes([0.2, 0.06, 0.6, 0.02])
        self.slider_x0 = Slider(ax_x0, 'Initial x₀ [m]',
                                0.01, 0.12,
                                valinit=self.osc.initial_position,
                                valstep=0.01)

        ax_reset = plt.axes([0.42, 0.01, 0.16, 0.03])
        self.btn_reset = Button(ax_reset, 'Reset Simulation')

        self.slider_mass.on_changed(self.on_parameter_change)
        self.slider_damping.on_changed(self.on_parameter_change)
        self.slider_x0.on_changed(self.on_parameter_change)
        self.btn_reset.on_clicked(self.on_reset)

    def on_parameter_change(self, val):
        self.osc.mass = self.slider_mass.val
        self.osc.d = self.slider_damping.val
        self.osc.initial_position = self.slider_x0.val
        self.update_info_display()

    def update_info_display(self):
        info = self.osc.compute_system_info()
        info_str = f"""m = {self.osc.mass:.2f} kg
k = {self.osc.k:.1f} N/m
d = {self.osc.d:.2f} Ns/m
ω₀ = {info['omega_0']:.2f}
D = {info['damping_ratio']:.3f}
{info['damping_type']}"""
        self.info_text.set_text(info_str)

    def on_reset(self, event):
        self.osc.reset()
        self.frame = 0
        self.position_history[:] = 0
        self.velocity_history[:] = 0

    def animate(self, frame_num):
        if self.frame >= self.steps:
            return

        x, v, a = self.osc.step(self.dt)

        self.position_history[self.frame] = x
        self.velocity_history[self.frame] = v

        mass_x = x - 0.0175
        self.mass_rect.set_xy((mass_x, 0.0825))

        spring_x, spring_y = self.draw_spring(-0.12, mass_x)
        self.spring_line.set_data(spring_x, spring_y)

        self.pos_line.set_data(self.time_array[:self.frame+1],
                               self.position_history[:self.frame+1])
        self.vel_line.set_data(self.time_array[:self.frame+1],
                               self.velocity_history[:self.frame+1])

        self.frame += 1
        return self.mass_rect, self.spring_line, self.pos_line, self.vel_line

    def run(self):
        self.ani = FuncAnimation(self.fig, self.animate,
                                 frames=self.steps,
                                 interval=self.dt * 1000,
                                 blit=False)
        plt.show()


# =====================================================
# Main
# =====================================================
if __name__ == "__main__":

    oscillator = SimpleOscillator(
        mass=0.5,
        spring_constant=20.0,
        damping=0.5,
        initial_position=0.08
    )

    visualizer = OscillatorVisualizer(oscillator)
    visualizer.run()