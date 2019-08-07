"""
Example of how to use the mirisim_tso package
"""

import os

import mirisim_tso
mirisim_tso.utils.init_log(stdout_loglevel="INFO", file_loglevel="DEBUG")

print("-------------------------------------------------------------------------------------------------------")
print("                                 POST-TREATMENT OF MIRISIM SIMULATIONS                                 ")
print("                                       TIME-SERIES OBSERVATIONS                                        ")
print("-------------------------------------------------------------------------------------------------------")

simulation_name = "MIRI_1Exposure2h"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
simulation_directory = os.path.join(simulation_name)
det_images_directory = os.path.join(simulation_directory, "det_images")
illum_models_directory = os.path.join(simulation_directory, "illum_models")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
t = 60  # seconds
x = 245  # pixel
y = 38  # pixel
frame = 0.19  # seconds

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
signal = mirisim_tso.utils.extract_signal(illum_models_directory, x, y)
original_ramp = mirisim_tso.utils.read_det_image(det_images_directory, x, y)
ramp_difference = mirisim_tso.effects.response_drift(original_ramp, t, signal, frame)
