import mirisim_tso
import numpy as np
import matplotlib.pyplot as plt


# import logging
# LOG = logging.getLogger(__name__)


def plot_response_drift():
    """
    Test response drift effect

    :return:
    """

    t_0 = 60  # seconds
    frame = 0.19  # seconds
    nb_frames = 25
    flux = 10.  # DN/s
    flux = np.array(flux)
    subarray_x_size = 72
    subarray_y_size = 416

    frame_sample = np.arange(nb_frames)
    original_ramp = frame_sample * flux * frame
    time_sample = frame_sample * frame

    original_ramp = original_ramp[np.newaxis,:,np.newaxis, np.newaxis]

    time = []
    time_flux = []
    for t_0 in np.arange(0, 12000, 100):
        ramp_difference = mirisim_tso.effects.response_drift(original_ramp, t_0, flux, frame)
        ramp_i = np.squeeze(original_ramp) + np.squeeze(ramp_difference)

        (fit_flux_i, dummy) = np.polyfit(time_sample, ramp_i, 1)

        time.append(t_0)
        time_flux.append(fit_flux_i)

    fig = plt.figure()
    ax = fig.add_subplot(2, 1, 1)
    ramp_difference = mirisim_tso.effects.response_drift(original_ramp, t_0, flux, frame)
    ax.plot(frame_sample, np.squeeze(original_ramp), label="Original Ramp")
    ax.plot(frame_sample, np.squeeze(original_ramp) + np.squeeze(ramp_difference), label="response drift applied")
    # ax.plot(frame_sample, ramp_difference, label="response drift difference")
    ax.legend()
    ax = fig.add_subplot(2, 1, 2)
    ax.plot([0, 12000], [flux, flux], label="Input flux")
    ax.plot(time, time_flux, label="RD corrected flux")
    ax.legend()

    fig.suptitle("Response drift")
    plt.savefig("response_drift.png")

    return fig


def plot_anneal_recovery():
    """
    Test response drift effect

    :return:
    """

    config = {'simulations': {'dir': 'MIRI_1Integration', 'overwrite': True},
              'response_drift': {'active': True},
              'idle_recovery': {'active': True, 'duration': 1000.},
              'anneal_recovery': {'active': True, 'time': -0.},
              'noise': {'active': True}}

    frame = 0.19  # seconds
    nb_frames = 25
    flux = 10  # DN/s

    frame_sample = np.arange(nb_frames)
    original_ramp = frame_sample * flux * frame
    time_sample = frame_sample * frame

    original_ramp = original_ramp[np.newaxis,:,np.newaxis, np.newaxis]

    time = []
    time_flux = []
    for t_0 in np.arange(0, 2000, 10):
        ramp_difference = mirisim_tso.effects.anneal_recovery(original_ramp, t_0, frame, config)
        ramp_i = np.squeeze(original_ramp) + np.squeeze(ramp_difference)

        (fit_flux_i, dummy) = np.polyfit(time_sample, ramp_i, 1)

        time.append(t_0)
        time_flux.append(fit_flux_i)

    fig = plt.figure()
    ax = fig.add_subplot(2, 1, 1)
    t_0 = 0
    ramp_difference = mirisim_tso.effects.anneal_recovery(original_ramp, t_0, frame, config)
    ax.plot(frame_sample, np.squeeze(original_ramp), label="Original Ramp")
    ax.plot(frame_sample, np.squeeze(original_ramp) + np.squeeze(ramp_difference), label="Anneal recovery applied with t_0={}".format(t_0))
    ax.set_xlabel("Frame number")
    ax.set_ylabel("Count [DN]")
    # ax.plot(frame_sample, ramp_difference, label="response drift difference")
    ax.legend()
    ax = fig.add_subplot(2, 1, 2)
    ax.plot([0, 2000], [flux, flux], label="Input flux")
    ax.plot(time, time_flux, label="AR corrected flux")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Flux [DN/s]")
    ax.legend()
    fig.suptitle("Anneal recovery")
    fig.savefig("anneal.png")

    return fig


def plot_idle_recovery():
    """
    Test response drift effect

    :return:
    """


    config = {'simulations': {'dir': 'MIRI_1Integration', 'overwrite': True},
              'response_drift': {'active': True},
              'idle_recovery': {'active': True, 'duration': 5000.},
              'anneal_recovery': {'active': True, 'time': -0.},
              'noise': {'active': True}}


    frame = 0.19  # seconds
    nb_frames = 25
    flux = 500  # DN/s

    frame_sample = np.arange(nb_frames)
    original_ramp = frame_sample * flux * frame
    time_sample = frame_sample * frame



    original_ramp = original_ramp[np.newaxis,:,np.newaxis, np.newaxis]

    time = []
    time_flux = []
    for t_0 in np.arange(0, 12000, 100):
        ramp_difference = mirisim_tso.effects.idle_recovery(original_ramp, t_0, flux, frame, config)
        ramp_i = np.squeeze(original_ramp) + np.squeeze(ramp_difference)

        (fit_flux_i, dummy) = np.polyfit(time_sample, ramp_i, 1)

        time.append(t_0)
        time_flux.append(fit_flux_i)

    fig = plt.figure()
    ax = fig.add_subplot(2, 1, 1)
    t_0 = 0.
    ramp_difference = mirisim_tso.effects.idle_recovery(original_ramp, t_0, flux, frame, config)
    ax.plot(frame_sample, np.squeeze(original_ramp), label="Original Ramp")
    ax.plot(frame_sample, np.squeeze(original_ramp) + np.squeeze(ramp_difference), label="Idle recovery applied, t_0={} s".format(t_0))
    # ax.plot(frame_sample, ramp_difference, label="response drift difference")
    ax.legend()

    ax = fig.add_subplot(2, 1, 2)
    ax.plot([0, 12000], [flux, flux], label="Input flux")
    ax.plot(time, time_flux, label="IR corrected flux")
    ax.legend()
    fig.suptitle("Idle recovery")
    plt.savefig("idle_recovery.png")

    return fig

if __name__ == "__main__":
    print("-------------------------------------------------------------------")
    print("|                           SCRIPT PLOT                           |")
    print("-------------------------------------------------------------------")
    #mirisim_tso.utils.init_log(stdout_loglevel="DEBUG", file_loglevel="DEBUG")
    #plot_anneal_recovery()
    plot_idle_recovery()

    plt.show()
