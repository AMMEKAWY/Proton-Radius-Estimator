import subprocess
import numpy as np

num_runs = 50
mean_results = np.zeros(num_runs)

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

        # Assign the proton radius to the corresponding position in the mean_results array
        mean_results[i] = round(float(output) * 1e15, 2)

        # Print the proton radius for the current run
        print(f"Proton radius of run {i + 1} = {mean_results[i]} fm")

overall_mean_radius = round(np.mean(mean_results), 2)

# Calculate the Standard Deviation
std_dev = np.std(mean_results, ddof=1)

# Print the mean proton radius and standard deviation
print(f"\nMean proton radius over {num_runs} runs = {overall_mean_radius} fm +/- {round(std_dev,2)} fm (1 SD)")

# Assume you have the true value for the proton radius
true_proton_radius = 0.84  # Replace with the actual true value
print("Actual proton radius =", true_proton_radius)

# Calculate the Relative error
re = abs(overall_mean_radius - true_proton_radius) / true_proton_radius
print(f"\nCalculated proton radius over {num_runs} runs = {overall_mean_radius} fm")
# Print the error results
print("Relative error =", round(re,2))
