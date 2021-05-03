import numpy as np

def aniso_parametrization(x, A, B, u_0, phi_2, phi_4):
    """
    Parametrization of azimuthal anisotropy, following Smith and Dahlen (1973).

    > Smith, M. L., & Dahlen, F. A. (1973). The azimuthal dependence of Love and Rayleigh wave propagation in a slightly anisotropic medium. Journal of Geophysical Research, 78(17), 3321–3333. http://doi.org/10.1029/JB078i017p03321

    :param x: Azimuth [rad]
    :param A: Peak-to-Peak Ampltiude of 2Theta-term
    :param B: Peak-to-Peak Ampltiude of 4Theta-term
    :param u_0: Mean Velocity / Absolute Shift
    :param phi_2: Fast axis orientation 2Theta-term [rad]
    :param phi_4: Fast axis orientation 4Theta-term [rad]
    """
    return u_0 + A * np.cos(2 * (x - phi_2)) + B * np.cos(4 * (x - phi_4))


def fit_aniso(az_bins, az_step, vel_bin_medians, vel_bin_stds):
    """ Fit 4-Theta anisotropy in azimuth-binned data. """
    from scipy import optimize
    import logging

    # initial guesses for anisotropic terms to stabilize curve_fit
    A_guess = np.random.uniform(
        low=np.mean(vel_bin_medians) + np.min(vel_bin_medians),
        high=np.mean(vel_bin_medians) + np.max(vel_bin_medians)
        )
    B_guess = np.random.uniform(
        low=np.mean(vel_bin_medians) + np.min(vel_bin_medians),
        high=np.mean(vel_bin_medians) + np.max(vel_bin_medians)
        )
    u0_guess = np.random.uniform(
        low=np.min(vel_bin_medians),
        high=np.max(vel_bin_medians)
        )
    phi2_guess = np.random.uniform(
        low=0,
        high=2*np.pi
        )
    phi4_guess = np.random.uniform(
        low=0,
        high=2*np.pi
        )
    p0 = (A_guess, B_guess, u0_guess, phi2_guess, phi4_guess)

    try:
        params, params_covariance = optimize.curve_fit(
            f=aniso_parametrization,
            xdata=np.deg2rad(az_bins + az_step/2),
            ydata=vel_bin_medians,
            sigma=vel_bin_stds,
            absolute_sigma=True,
            p0=p0)
    except RuntimeError:
        logging.warning('skipping, fit failed')

    return (params, params_covariance)


def extract_aniso_results(params, params_covariance):
    param_std = np.sqrt(np.diag(np.abs(params_covariance)))

    # extract results from anisotropic fit
    A, B, u_0, phi_2, phi_4 = params
    A_std, B_std, u_0_std, phi_2_std, phi_4_std = param_std

    # relative amplitudes and errors
    rel_A = 100 * A / u_0
    rel_A_std = 100 * A_std / rel_A
    rel_B = 100 * B / u_0
    rel_B_std = 100 * B_std / rel_B

    # correct phi_2 if negative amplitude A was fitted
    phi_2 = np.degrees(phi_2)
    if A < 0:
        phi_2 += 90
    phi_2 %= 180

    phi_2_std = np.degrees(phi_2_std)

    return phi_2, phi_2_std, rel_A, rel_A_std
