import mirisim_tso
import numpy as np
import pytest
import numpy.testing as npt

testdata = [  # flux in DN/s, t_0 in s, frame_time time in s, mean ramp (DN)
    (4999., 10000., 0.19,  np.array([2, 3, 4, 5])),
]

@pytest.mark.parametrize("flux, t_0, frame_time, input_ramp", testdata)
def test_response_drift(input_ramp, flux, t_0, frame_time):
    """
    Test response drift effect

    We ensure that response drift converges, hence doesn't contribute past a certain time (here 10000s) at a certain
    flux (here 4999 DN/s).

    Parameters
    ----------
    input_ramp: np.ndarray
        ramp from ground testing
    flux: float
        mean flux in DN/s associated with the input
    t_0: float
        start time of the ramp, compared to the start of the observation (s)

    Returns
    -------

    """

    nb_frames = input_ramp.size

    frame_sample = np.arange(nb_frames)
    original_ramp = frame_sample * flux * frame_time

    transformed_ramp = original_ramp[np.newaxis, :, np.newaxis, np.newaxis]
    transformed_flux = np.full((1,1), flux)

    ramp_difference = mirisim_tso.effects.response_drift(transformed_ramp, t_0, transformed_flux, frame_time)

    # Note that the reference ramp is the difference between the original ramp, and its first element (offset substraction)
    npt.assert_array_almost_equal(ramp_difference, np.zeros_like(ramp_difference))


