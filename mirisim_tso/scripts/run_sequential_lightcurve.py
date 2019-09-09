#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Example of how to use the mirisim_tso package
"""

import mirisim_tso

lightcurve_folder = "../../test_light_curve"

config_filename = "../../post_treatment.ini"

mirisim_tso.sequential_lightcurve_post_treatment(lightcurve_folder, conf=config_filename)
