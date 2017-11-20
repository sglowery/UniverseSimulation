import matplotlib.pyplot as plt
from universe import Universe
import numpy as np
import conversions
import constants

benchmark = Universe(.31, .69, 9e-5, 68)

print(benchmark.density_at_scale(1))

print("Current age of the Benchmark model (time at scale factor = 1):",
      np.round(benchmark.cosmic_time(1, 1e3),decimals=2), "Gyr")
