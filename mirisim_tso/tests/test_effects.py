import mirisim_tso
import numpy as np

def test_response_drift():
    """
    Test response drift effect

    :return:
    """

    import matplotlib.pyplot as plt

    t_0 = 3600  # seconds
    x = 245  # pixel
    y = 38  # pixel
    frame = 0.19  # seconds


    nb_frames = 25
    flux = 960  # DN/s
    frame_sample = np.arange(nb_frames)
    original_ramp = frame_sample * flux
    time_sample = frame_sample * frame

    time = []
    time_flux = []
    for t_0 in np.arange(0, 5400, 500):

        ramp_difference = mirisim_tso.effects.response_drift(original_ramp, t_0, flux, frame)
        ramp_i = original_ramp + ramp_difference

        (fit_flux_i, dummy) = np.polyfit(time_sample, ramp_i, 1)

        time.append(t_0)
        time_flux.append(fit_flux_i)


        #time.extend(time_sample)
        #ramps.extend(ramp_i)


    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    #ax.plot(frame_sample, original_ramp, label="Original Ramp")
    #ax.plot(frame_sample, original_ramp+ramp_difference, label="response drift applied")
    #ax.plot(frame_sample, ramp_difference, label="response drift difference")
    ax.plot(time, time_flux, label="Flux over time")
    ax.plot([0, 5400], [flux, flux], label="Input flux")

    ax.legend()
    plt.show()