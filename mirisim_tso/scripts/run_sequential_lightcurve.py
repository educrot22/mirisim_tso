#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Example of how to use the mirisim_tso package
"""

import mirisim_tso

config_filename = "../../post_treatment.ini"

mirisim_tso.sequential_lightcurve_post_treatment(conf=config_filename)
