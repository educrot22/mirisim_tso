#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Example of how to use the mirisim_tso package
"""

import mirisim_tso

mirisim_tso.utils.init_log(stdout_loglevel="INFO", file_loglevel="DEBUG")

simulation_name = "./simulation_870"

config_filename = "post_treatment.ini"

mirisim_tso.single_simulation_post_treatment(simulation_folder=simulation_name, t_0=0., conf=config_filename)
