#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Example of how to use the mirisim_tso package
"""
import sys
import mirisim_tso

if len(sys.argv) > 1:
    config_filename = sys.argv[1]
else:
    sys.stdout.write('Configuration file missing')
    sys.exit()

mirisim_tso.sequential_lightcurve_post_treatment(config_filename)

