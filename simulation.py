import matplotlib.pyplot as plt
from universe import Universe
import numpy as np

# First I create an array and use numpy to populate the list with times, starting at 1 second and going to 500 Gyr.
# Obviously, the array's steps are logarithmic.
time_list = np.logspace(np.log10(1), np.log10(500*1e9*365*24*3600), num=1e5)
# Empty lists to contain values of component densities at different scale factors
matter_densities_list = []
lambda_densities_list = []
radiation_densities_list = []
scales_list = []

benchmark = Universe(.31, 9e-5, .69, 68)
# First, calculate how scale factor changes over time and put each entry into the scales_list array
for time in time_list:
    scales_list.append(benchmark.scale_at_time(time))

# This creates the scale vs. time plot with a vertical, black, dashed line at 13.75 Gyr
plt.plot(time_list, scales_list, label="Scale factor")
plt.xlabel("Time (sec)")
plt.ylabel("Scale factor")
plt.xscale("log")
plt.yscale("log")
plt.title("Scale factor vs. Time")
plt.axvline(x=13.75*1e9*365*24*3600, color='black', linestyle='--')
plt.show()

# Now, calculate component densities using the previously generated list of scale factors
for scale in scales_list:
    matter_density, lambda_density, radiation_density = benchmark.density_at_scale(scale)

    matter_densities_list.append(matter_density)
    lambda_densities_list.append(lambda_density)
    radiation_densities_list.append(radiation_density)

# Plot the three components and another vertical, black, dashed line at a=1
plt.plot(scales_list, matter_densities_list, label="Matter")
plt.plot(scales_list, lambda_densities_list, label="Cosmological constant")
plt.plot(scales_list, radiation_densities_list, label="Radiation")
plt.axvline(x=1, color='black', linestyle='--')
plt.legend(loc="center left", shadow=True, fontsize="small")
plt.xlabel("Scale factor")
plt.ylabel("Density parameter")
plt.xscale("log")
plt.title("Density vs. scale factor")
plt.show()
