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

simulation_name = "MIRI_1Integration"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
simulation_directory = os.path.join(simulation_name)
det_images_filename = os.path.join(simulation_directory, "det_images", "det_image_seq1_MIRIMAGE_P750Lexp1.fits")
illum_models_filename = os.path.join(simulation_directory, "illum_models", "illum_model_1_MIRIMAGE_P750L.fits")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
t_0 = 60  # seconds
frame_time = 0.19  # seconds

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

signal = mirisim_tso.utils.read_illum_model(illum_models_filename)
original_ramp = mirisim_tso.utils.read_det_image(det_images_filename)
ramp_difference = mirisim_tso.effects.response_drift(original_ramp, t_0, signal, frame_time)

new_ramp = original_ramp + ramp_difference

metadatas = {'history':"Post processing with MIRISim TSO v{}".format(mirisim_tso.__version__)}



# Write fits file
mirisim_tso.utils.write_det_image_with_effects(det_images_filename, new_data=new_ramp, extra_metadata=metadatas)




