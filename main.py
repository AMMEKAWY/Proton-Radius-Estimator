import subprocess
import numpy as np

num_runs = 10
mean_results = []

for i in range(num_runs):
    # Run the sps.py script using subprocess
    process = subprocess.Popen(['python3', 'sps.py'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    output, errors = process.communicate()

    # Check if the process completed successfully
    if process.returncode == 0:

        # Extract the proton radius from the output (modify this based on the actual output format)
        
        lines = output.splitlines()
        for line in lines:
            if "Proton radius of run" in line:
                try:
                    proton_radius = float(line.split('=')[1].strip())
                    break
                except ValueError:
                    print(f"Error: Invalid value for proton radius in output for run {i + 1}")

            mean_results.append(float(output) * 1e15)

            # Print the proton radius for the current run
            print(f"Proton radius of run {i + 1} = {float(output) * 1e15}")


overall_mean_radius = np.mean(mean_results)
# Print the mean proton radius
print(f"\nThe mean proton radius over {num_runs} runs = {overall_mean_radius}")

# Assume you have the true value for the proton radius
true_proton_radius = 0.84  # Replace with the actual true value

# Calculate the absolute error
absolute_error = abs(overall_mean_radius - true_proton_radius)

# Calculate the relative error as a percentage
relative_error_percentage = (absolute_error / true_proton_radius) * 100

# Print the error results
print(f"Measured Proton Radius: {overall_mean_radius}")
print(f"True Proton Radius: {true_proton_radius}")
print(f"Absolute Error: {absolute_error}")
print(f"Relative Error (%): {relative_error_percentage}")

