import matplotlib.pyplot as plt
from universe import Universe
import numpy as np
import conversions
import constants

benchmark = Universe(.31, .69, 9e-5, 68)

print(benchmark.density_at_scale(10**-5))
