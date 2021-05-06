#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Example of how to use the mirisim_tso package
"""
import time
import sys
#
import mirisim_tso
#
if len(sys.argv) > 1:
    config_filename= sys.argv[1]
else:
    print("Error: You need to provide an .ini filename as argument of the script")
    sys.exit()
#
start_time = time.time()
print('mirisim_tso', mirisim_tso.__version__, config_filename)
mirisim_tso.utils.init_log(stdout_loglevel="INFO", file_loglevel="DEBUG")
mirisim_tso.sequential_lightcurve_post_treatment(config_filename)
elapsed_time = time.time() - start_time
print('elapsed time:', elapsed_time)
