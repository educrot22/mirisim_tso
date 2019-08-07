import mirisim_tso
import numpy as np
import pytest
import numpy.testing as npt

testdata = [  # flux in DN/s, t_0 in s, frame time in s, mean ramp (DN)
    (500., 3600., 0.19,  np.array([2, 3, 4, 5])),
]

@pytest.mark.parametrize("flux, t_0, frame_time, input_ramp", testdata)
def test_response_drift(input_ramp, flux, t_0, frame_time):
    """
    Test response drift effect

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

    ramp_difference = mirisim_tso.effects.response_drift(original_ramp, t_0, flux, frame_time)
    ramp_i = original_ramp + ramp_difference

    # We can't compare the first frame because of the ramp offset in test data.
    # So first frame will always be perfectly 0 in test data
    ref_ramp = input_ramp[1:] - input_ramp[0]

    ramp_i = ramp_i[1:]

    npt.assert_approx_equal(ramp_i, ref_ramp)


