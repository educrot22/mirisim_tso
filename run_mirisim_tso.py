"""
Example of how to use the mirisim_tso package
"""

import os

import mirisim_tso

simulation_name = "MIRI_1Integration"

config_filename = "post_treatment.ini"

mirisim_tso.single_simulation_post_treatment(simulation_folder=simulation_name, t_0=1000., conf=config_filename)
