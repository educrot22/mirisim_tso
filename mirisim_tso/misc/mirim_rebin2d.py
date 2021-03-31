from __future__ import print_function
import numpy as np
import logging

LOG = logging.getLogger('kiss.mirim_rebin2d')


def mirim_rebin2d(data, zoom):
    """
    rebin as idl
    
    TODO RG
    
    Parameters
    ~~~~~~~~~~
    :param nd.array data: numpy array of any dimensions shape
    
    :param float zoom: ?

    :return: numpy array shrinked to shape
    
    Example
    ~~~~~~~
    >>> a = np.arange(48).reshape((6,8))
    >>> mirim_rebin2d( a, [2,2])
        
        
    HISTORY:
    -------------
    Rene Gastaud, 13 January 2016
    https://codedump.io/share/P3pB13TPwDI3/1/resize-with-averaging-or-rebin-a-numpy-2d-array
    2017-11-17 : RG and Alan O'Brien  compatibility with python 3 bug452
    """
    #  rebin factor
    # sh = data.shape[0]/zoom[0], zoom[0], data.shape[1]/zoom[1], zoom[1]
    #  just in case  Belt and suspenders
    sh = int(data.shape[0] // zoom[0]), zoom[0], int(data.shape[1] // zoom[1]), zoom[1]
    LOG.debug("shape: {}".format(sh))
    result = data.reshape(sh).sum(3).sum(1)
    return result
