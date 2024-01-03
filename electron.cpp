#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include <unordered_map>
#include <algorithm>

int n=500000;

// Constants
const double q_e = -1.602176634e-19;   // Elementary charge (C)
const double epsilon_0 = 8.854187817e-12;  // Vacuum permittivity (F/m)
const double m_e = 9.10938356e-31;      // Electron mass (kg)
const double h_bar = 1.054571817e-34;   // Reduced Planck constant (Js)
std::vector<double> angle;
//const std::string filename = "deflection_angles.txt";
//std::ofstream file(filename, std::ios::app);

// Function to calculate Coulomb force in x and y directions
void coulombForce(double q1, double q2, double x, double y, double z, double& fx, double& fy, double& fz) {
    double r = sqrt(x * x + y * y + z * z);
    double coeff = (q1 * q2) / (epsilon_0 * r * r);

    
    fx = coeff * x / r;  // Coulomb force in x direction
    fy = coeff * y / r;  // Coulomb force in y direction
    fz = coeff * z / r;  // Coulomb force in z direction

}

// Function to update particle position and velocity using the Euler method
void updateParticle(double& x, double& y, double& z, double& vx, double& vy, double& vz, double fx, double fy, double fz, double dt) {
    // Update velocity
    vx += fx * dt / m_e;
    vy += fy * dt / m_e;
    vz += fz * dt / m_e;

    //std::cout<<vx<<std::endl;

    // Update position
    x += vx * dt;
    y += vy * dt;
    z += vz * dt;
}
/*
// Function to write data to a file
void writeDataToFile(const std::string& filename, double t) {

    file << t <<"\n";
    
}*/

// Function to perform the simulation with specified initial conditions
void simulate(double mean_position, double delta_y, double mean_momentum, double dt, double t_max) {

    std::random_device rd;
    std::default_random_engine generator(rd());

    // Set the desired mean and standard deviation
    double mean = mean_position;
    double stddev = delta_y;

    // Create a normal distribution with the specified mean and standard deviation
    std::normal_distribution<double> distribution(mean, stddev);
    std::uniform_real_distribution<double> momentum_distribution(0, 2*mean_momentum);


    // Sample initial conditions
    double initial_x = -1e-7;
    double initial_y = distribution(generator);
    double initial_z = 0.0;
    
    double momentum_sample = momentum_distribution(generator);
    double initial_vx = momentum_sample / m_e;
    double initial_vy = 0.0;
    double initial_vz = 0.0;
    
    double last_x[4], last_y[4];
    
    // Simulation loop
    for (double t = 0; t <= t_max; t += dt) {
        // Calculate Coulomb forces in x, y, and z directions
        double fx, fy, fz;
        coulombForce(q_e, q_e, initial_x, initial_y, initial_z, fx, fy, fz);

        // Update particle position and velocity
        updateParticle(initial_x, initial_y, initial_z, initial_vx, initial_vy, initial_vz, fx, fy, fz, dt);

        double progress = t / t_max;

        // Update the arrays with the last four x, y values
        last_x[0] = last_x[1];
        last_x[1] = last_x[2];
        last_x[2] = last_x[3];
        last_x[3] = initial_x;

        last_y[0] = last_y[1];
        last_y[1] = last_y[2];
        last_y[2] = last_y[3];
        last_y[3] = initial_y;
                    
    }
    
    
    double dy=last_y[3]-last_y[2];
    double dx=last_x[3]-last_x[2];
    
    double angle_value= std::atan2(dy, dx);
    
    if(angle_value<0){ angle_value = (angle_value)+2*M_PI;}
    
    angle.push_back(angle_value);
    //std::cout<<angle<<std::endl;
    //writeDataToFile(filename,angle_value);
}



void createHistogram(const std::vector<double>& data, int numBins) {

    std::ofstream file2("H_50.txt");

    double minRange = 0.0;
    double maxRange = M_PI;  // Using M_PI constant from cmath for pi
    double binWidth = (maxRange - minRange) / numBins;

    for (int i = 0; i < numBins; ++i) {
        double binStart = minRange + i * binWidth;
        double binEnd = binStart + binWidth;

        int count = 0;

        // Count the elements within the current bin
        for (double num : data) {
            if (num >= binStart && num < binEnd) {
                count++;
            }
        }

        // Calculate kappa as a ratio of the count to the total number of elements
        double kappa = static_cast<double>(count) / data.size();

        // Display asterisks based on the count
        file2 << (binStart+binEnd)/2 <<" "<< kappa << std::endl;
        //std::cout << (binStart+binEnd)/2 <<" "<< kappa << std::endl;
    }
    
    file2.close();
}


int main() {

    

    // Parameters based on specified formulas
    double mean_energy = 50 * std::abs(q_e);
    double mean_momentum = std::sqrt(2 * m_e * mean_energy);
    double mean_position = 5.29e-11;
    
    std::random_device rd;
    std::default_random_engine generator(rd());

    std::uniform_real_distribution<double> y_distribution(-mean_position, mean_position);

    double delta_y = h_bar / (2 * std::abs(mean_momentum / std::sqrt(12)));
    double dt = 1.0e-17;  // Time step (s)
    double t_max = 1.0e-13;  // Maximum simulation time (s)
    
    for(int i=0; i<n; i++){
    
    double randy=y_distribution(generator);
    
    // Call the simulate function with specified initial conditions
    simulate(mean_position, randy, mean_momentum, dt, t_max);

    double j=(float)(i + 1) / n * 100.0;
    
    std::cout << "Simulation progress: " << j << "%" << std::flush;
    std::cout << "\r";  // Move the cursor back to the beginning of the line
 

    }
    
    std::cout<<"			"<<std::endl;
    //file.close();
    std::cout<<"creating Histogram"<<std::endl;
    createHistogram(angle, 180);

    
    std::cout<<"done!"<<std::endl;

    return 0;
}

