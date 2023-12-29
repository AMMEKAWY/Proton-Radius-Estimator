import numpy as np
import matplotlib.pyplot as plt

angles = []
deflection_points = []

def coulomb_force(q1, q2, r):
    k = 8.9875e9  # Coulomb's constant in N m^2 / C^2
    return k * (q1 * q2) / r**2

def simulate_motion(charge, mass, initial_position, initial_velocity, num_steps, dt):
    positions = np.zeros((num_steps, 2))
    velocities = np.zeros((num_steps, 2))
    angles = []

    positions[0] = initial_position
    velocities[0] = initial_velocity



    for step in range(1, num_steps):
        x, y = positions[step - 1]
        vx, vy = velocities[step - 1]

        r = np.sqrt(x**2 + y**2)

        # Calculate Coulomb force components
        fx = coulomb_force(charge, -1.0, r) * x / r
        fy = coulomb_force(charge, -1.0, r) * y / r

        # Update velocity using F = ma
        ax = fx / mass
        ay = fy / mass
        vx += ax * dt
        vy += ay * dt

        # Update position using velocity
        x += vx * dt
        y += vy * dt

        positions[step] = [x, y]
        velocities[step] = [vx, vy]

        angle = np.arctan2(y - positions[step - 1, 1], x - positions[step - 1, 0])


        if angle < 0:
            angle += np.pi * 2
        
            angles.append(angle)

        if round(angle, 1) > np.pi/2 or round(angle, 1) < 3*np.pi/2:
            deflection_points.append(positions[0, 1])

    return positions, angles, deflection_points

def run_simulation_samples(num_samples_y, num_samples_p, charge, mass, mean_position, mean_momentum, num_steps, dt):
    for i in range(num_samples_y):
        # Generate a random position from a normal distribution based on Heisenberg uncertainty
        delta_y = h_bar / (2 * np.abs(mean_momentum / np.sqrt(12)))
        sampled_position = np.random.normal(loc=mean_position, scale=delta_y)
        
        for j in range(num_samples_p):
            # Generate a random momentum from a uniform distribution
            momentum_sample = np.random.uniform(0, 6e6 * mean_momentum)

            # Run simulation for the current position and momentum sample
            initial_velocity_sample = momentum_sample / mass
            simulate_and_plot(charge, mass, [0.0, sampled_position], [initial_velocity_sample, 0.0], num_steps, dt)

def simulate_and_plot(charge, mass, initial_position, initial_velocity, num_steps, dt):
    particle_positions, angles, deflection_points = simulate_motion(charge, mass, initial_position, initial_velocity, num_steps, dt)
    plot_motion(particle_positions)

def calculate_differential_cross_section(scattering_angles, bin_width):
    # Create a histogram
    hist, bin_edges = np.histogram(scattering_angles, bins=np.arange(0, np.pi + bin_width, bin_width))

    # Calculate differential cross section for each bin
    total_particles = len(scattering_angles)
    delta_omega = bin_width
    differential_cross_section = hist / (total_particles * delta_omega)

    print(np.mean(angles) * 180 / np.pi, "		", np.max(differential_cross_section) / (2 * np.pi))

    return bin_edges, differential_cross_section

def plot_motion(positions):
    
    plt.plot(positions[:, 0], positions[:, 1])
    plt.title('Motion of Electron around Hydrogen Atom')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    

# Set parameters for an electron and a hydrogen atom
h_bar = 1.04e-34
energy = 2e7
electron_charge = -1.602e-19  # Charge of the electron in C
electron_mass = 9.109e-31  # Mass of the electron in kg
mean_momentum = np.sqrt(2 * electron_mass * energy * abs(electron_charge))  # Mean momentum in kg m/s
mean_position = 1e-19
initial_electron_position = [0.0, mean_position]  # Initial position of the electron in meters

num_steps = 50
dt = 1e-14  # Time step in seconds
num_samples_y = 50
num_samples_p = num_samples_y

# Run simulations for different sampled positions and momenta
run_simulation_samples(num_samples_y, num_samples_p, electron_charge**2, electron_mass, mean_position, mean_momentum, num_steps, dt)

def calculate_differential_cross_section(scattering_angles, bin_width):
    # Create a histogram
    hist, bin_edges = np.histogram(scattering_angles, bins=np.arange(min(scattering_angles), max(scattering_angles) + bin_width, bin_width))
    
    # Find the bins corresponding to back deflection
    deflection_bins = np.where(np.abs(bin_edges) < np.pi/2)[0]

    # Sum up the histogram values for the deflection bins
    deflection_count = np.sum(hist[deflection_bins])
    
    c=(deflection_count/num_samples_y)
    
    # Calculate differential cross section for each bin
    total_particles = len(scattering_angles)
    delta_theta = bin_width
    differential_cross_section = hist * 2*np.pi / (total_particles * delta_theta)
    
    print(np.abs(np.mean(deflection_points)/c))
    
    # Calculate bin centers
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    return bin_centers, differential_cross_section

def plot_histogram_and_cross_section(bin_centers, differential_cross_section):
    # Plot the histogram and differential cross section
    plt.figure(figsize=(10, 4))

    plt.subplot(1, 2, 1)
    plt.hist(scattering_angles, bins=np.arange(min(scattering_angles), max(scattering_angles), bin_width), alpha=0.7)
    plt.title('Scattering Angle Histogram')
    plt.xlabel('Scattering Angle (radians)')
    plt.ylabel('Number of Events')

    plt.subplot(1, 2, 2)
    plt.plot(bin_centers, differential_cross_section, marker='o')
    plt.title('Differential Cross Section')
    plt.xlabel('Scattering Angle (radians)')
    plt.ylabel(r'$\frac{d\sigma}{d\theta}$')

    plt.tight_layout()
    #plt.show()


# Example data (replace this with your actual scattering angles)
scattering_angles = np.random.uniform(0, np.pi, size=1000)

# Choose bin width
bin_width = 0.1

# Calculate and plot the histogram and differential cross section
bin_centers, differential_cross_section = calculate_differential_cross_section(scattering_angles, bin_width)
plot_histogram_and_cross_section(bin_centers, differential_cross_section)

