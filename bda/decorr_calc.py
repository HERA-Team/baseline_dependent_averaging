# -*- coding: utf-8 -*-
# Copyright (c) 2018 Paul La Plante
# Licensed under the 2-clause BSD License
"""Supporting tools for BDA calculation."""

import numpy as np
from astropy import constants as const
from astropy.coordinates import Angle
from astropy import units

# define some HERA-specific constants
hera_latitude = Angle("-30:43:17.5", units.deg)


def _dudt(lx, ly, hour_angle, earth_omega, wavelength):
    """Compute the change in u due to the Earth's rotation.

    This function computes the change in the u-coordinate of a baseline due to
    the rotation of the Earth. Taken from Equation (42) in Wijnholds et al.
    (2018).

    Parameters
    ----------
    lx : astropy Quantity
        The length of the baseline along the x-axis in ECEF coordinates. Must be
        compatible with units of length.
    ly : astropy Quantity
        The length of the baseline along the y-axis in ECEF coordinates. Must be
        compatible with units of length.
    hour_angle : astropy Angle
        The hour angle at which to compute decorrelation. hour_angle == 0
        corresponds to zenith.
    earth_omega : astropy Quantity
        The Earth's rotation speed. Must be compatible with units of angle/time.
    wavelength : astropy Quantity
        The wavelength of the radiation in question. Must be compatible with
        units of length.

    Returns
    -------
    astropy Quantity
        The change in u in a quantity with units of angle/time.
    """
    hour_angle = max(
        min(hour_angle, Angle(90.0, units.degree)), Angle(-90.0, units.degree)
    )
    return (
        (lx * np.cos(hour_angle) - ly * np.sin(hour_angle)) * earth_omega / wavelength
    )


def _dvdt(lx, ly, hour_angle, dec, earth_omega, wavelength):
    """Compute the change in v due to the Earth's rotation.

    This function computes the change in the v-coordinate of a baseline due to
    the rotation of the Earth. Taken from Equation (42) in Wijnholds et al.
    (2018).

    Parameters
    ----------
    lx : astropy Quantity
        The length of the baseline along the x-axis in ECEF coordinates. Must be
        compatible with units of length.
    ly : astropy Quantity
        The length of the baseline along the y-axis in ECEF coordinates. Must be
        compatible with units of length.
    hour_angle : astropy Angle
        The hour angle at which to compute decorrelation. hour_angle == 0
        corresponds to zenith.
    dec : astropy Angle
        The declination at which to compute decorrelation. Note that declination
        depends on telescope location, so a source observed at local zenith
        would have a declination equal to that of the telescope.
    earth_omega : astropy Quantity
        The Earth's rotation speed. Must be compatible with units of angle/time.
    wavelength : astropy Quantity
        The wavelength of the radiation in question. Must be compatible with
        units of length.

    Returns
    -------
    astropy Quantity
        The change in u in a quantity with units of angle/time.
    """
    hour_angle = max(
        min(hour_angle, Angle(90.0, units.degree)), Angle(-90.0, units.degree)
    )
    dec = max(min(dec, Angle(90.0, units.degree)), Angle(-90.0, units.degree))
    return (
        (lx * np.sin(dec) * np.sin(hour_angle) + ly * np.sin(dec) * np.cos(hour_angle))
        * earth_omega
        / wavelength
    )


def decorr_pre_fs_int_time(frequency, baseline, pre_fs_int_time):
    """Compute the decorrelation due to integrating before fringe stopping.

    This function computes the amount of decorrelation that occurs before fringe
    stopping. This decorrelation is unavoidable and is due to the motion of the
    baseline in the uvw-plane before fringe stopping/rephasing to the phase
    center of the observation.

    Parameters
    ----------
    frequency : astropy Quantity
        The frequency of the radiation in question. Must be compatible with
        units of inverse time.
    baseline : astropy Quantity
        The length of the baseline. Must be compatible with units of length.
    pre_fs_int_time : astropy Quantity
        The amount of time integrated prior to fringe stopping. Must be
        compatible with units of time.

    Returns
    -------
    float
        The decorrelation fraction for the baseline, wavelength, and integration
        time specified.
    """
    wavelength = const.c / frequency.to(1 / units.s)
    earth_rot_speed = (Angle(360, units.deg) / units.sday).to(units.arcminute / units.s)
    max_resolution = Angle(np.arcsin(wavelength / baseline), units.rad)
    return float(pre_fs_int_time * earth_rot_speed / max_resolution.to(units.arcminute))


def decorr_chan_width(chan_width, baseline, corr_fov):
    """Compute the decorrelation due to the channel width.

    Parameters
    ----------
    chan_width : astropy Quantity
        The channel width. Must be compatible with units of inverse time.
    baseline : astropy Quantity
        The baseline length. Must be compatible with units of length.
    corr_fov : astropy Angle
        The opening angle of interest.

    Returns
    -------
    float
        The decorrelation fraction for the channel width, baseline, and opening
        angle specified.
    """
    return float(
        chan_width.to(1 / units.s) * baseline * np.sin(corr_fov.to(units.rad)) / const.c
    )


def decorr_post_fs_int_time(
    lx, ly, post_fs_int_time, corr_fov, frequency, telescope_latitude=hera_latitude
):
    """Compute the decorrelation from averaging together spectra from different times.

    This is computed by implementing Equation 41 from Wijnholds et al. (2018).

    Parameters
    ----------
    lx : astropy Quantity
        The length of the baseline along the x-axis in ECEF coordinates. Must be
        compatible with units of length.
    ly : astropy Quantity
        The length of the baseline along the y-axis in ECEF coordinates. Must be
        compatible with units of length.
    post_fs_int_time : astropy Quantity
        The length of post fringe-stopped integration for the baseline. This is
        the "fundamental time" that is used for computing the decorrelation.
        Must be compatible with units of time.
    corr_fov : astropy Angle
        The opening angle at which the maximum decorrelation is to be
        calculated. Because a priori it is not known in which direction the
        decorrelation will be largest, the expected decorrelation is computed in
        all 4 cardinal directions at `corr_fov_angle` degrees off of zenith,
        and the largest one is used. This is a "worst case scenario"
        decorrelation.
    frequency : astropy Quantity
        The frequency for which to compute the decorrelation. Must be compatible
        with units of inverse time.
    telescope_latitude : astropy Angle, optional
        The latitude for the telescope. This is needed to compute the change of
        the v-coordinate with respect to time due to the Earth's rotation. If
        not provided, the default is to use that of HERA.

    Returns
    -------
    decorr_frac : float
        The resulting decorrelation fraction for the specified parameters.
    rfac : float
        The amplitude reduction factor owing to motion in the uv-plane.
    """
    wavelength = const.c / frequency.to(1 / units.s)
    earth_rot_speed = (Angle(360, units.deg) / units.sday).to(units.arcminute / units.s)

    # case 1: +l
    du = _dudt(lx, ly, corr_fov, earth_rot_speed, wavelength)
    lval = np.cos(90.0 * units.deg + corr_fov)
    rfac = (du * lval) ** 2

    # case 2: -l
    du = _dudt(lx, ly, -corr_fov, earth_rot_speed, wavelength)
    lval = np.cos(90.0 * units.deg - corr_fov)
    rfac = max(rfac, (du * lval) ** 2)

    # case 3: +m
    dv = _dvdt(lx, ly, 0.0, telescope_latitude + corr_fov, earth_rot_speed, wavelength)
    mval = np.cos(90.0 * units.deg + corr_fov)
    rfac = max(rfac, (dv * mval) ** 2)

    # case 4: -m
    dv = _dvdt(lx, ly, 0.0, telescope_latitude - corr_fov, earth_rot_speed, wavelength)
    mval = np.cos(90.0 * units.deg - corr_fov)
    rfac = max(rfac, (dv * mval) ** 2)

    # make sure we have the right units
    rfac = rfac.to(units.rad ** 2 / units.s ** 2)

    # add other factors; return decorrelation fraction and max rfac value
    decorr_frac = (
        np.pi ** 2 * (post_fs_int_time.to(units.s).value) ** 2 / 6.0 * rfac.value
    )
    return decorr_frac, rfac.value


def bda_compression_factor(
    max_decorr=0.1,
    frequency=(250.0 * 1e6 * units.Hz),
    lx=(14.6 * units.m),
    ly=(14.6 * units.m),
    corr_fov_angle=Angle(20.0, units.degree),
    chan_width=(30.517 * units.kHz),
    pre_fs_int_time=(0.1 * units.s),
    corr_int_time=(10 * units.s),
):
    """Compute the total BDA compression factor for a given observation.

    This function finds the maximum number of power-of-two foldings that can be
    accommodated while staying below a specified maximum amount of decorrelation
    for an observation, taking into account decorrelation due to pre-fringe-
    stopping averaging, channel-width incoherence, and post-fringe-stopping
    averaging. This is useful for figuring out how many observations to combine
    to simulate the effects of BDA.

    Parameters
    ----------
    max_decorr : float
        The maximum amount of decorrelation desired. The calculated number of
        samples that can be averaged will be less than the specified level.
    frequency : astropy Quantity, optional
        The frequency at which to compute the decorrelation. Default is 250 MHz.
        If specified, must be compatible with units of inverse time.
    lx : astropy Quantity, optional
        The length of the baseline along the x-axis in ECEF coordinates. If
        specified, must be compatible with units of length.
    ly : astropy Quantity, optional
        The length of the baseline along the y-axis in ECEF coordinates. If
        specified, must be compatible with units of length.
    corr_fov_angle : astropy Angle, optional
        The opening angle at which the maximum decorrelation is to be
        calculated. Because a priori it is not known in which direction the
        decorrelation will be largest, the expected decorrelation is computed in
        all 4 cardinal directions at `corr_fov_angle` degrees off of zenith,
        and the largest one is used. This is a "worst case scenario"
        decorrelation.
    chan_width : astropy Quantity, optional
        The channel width to use for decorrelation owing to channel-width
        incoherence. If specified, must be compatible with units of inverse-
        time.
    pre_fs_int_time : astropy Quantity, optional
        The amount of time integrated prior to fringe stopping. If specified,
        must be compatible with units of time.
    corr_int_time : astropy Quantity, optional
        The length of post fringe-stopped integration for the baseline. This is
        the "fundamental time" that is used for computing the decorrelation.
        If specified, must be compatible with units of time.

    Returns
    -------
    num_two_foldings : int
        The number of power-of-two foldings that can be accommodated while
        remaining under the specified level of decorrelation. Should be
        interpreted as the power to which to raise 2 to get the total
        integration time. For example, if no additional averaging can be done,
        the value returned is 0. 2**0 == 1, and so only 1 time sample should be
        averaged together. If the value returned is 1, then 2**1 == 2 samples
        can be averaged, etc.
    """
    # calculate the pre-BDA decorrelation given the correlator settings
    baseline = np.sqrt(lx ** 2 + ly ** 2)
    decorr_cw = decorr_chan_width(chan_width, baseline, corr_fov_angle)

    decorr_pre_int = decorr_pre_fs_int_time(frequency, baseline, pre_fs_int_time)

    decorr_post_int, max_rfac = decorr_post_fs_int_time(
        lx, ly, corr_int_time, corr_fov_angle, frequency
    )

    # calculate pre- and post-fs decorrelations
    pre_fs_decorr = 1 - (1 - decorr_cw) * (1 - decorr_pre_int)
    total_decorr = 1 - (1 - pre_fs_decorr) * (1 - decorr_post_int)

    if total_decorr < max_decorr:
        # If total decorrelation is less than max allowed, we can average.
        # Figure out the maximum amount of decorrelation allowed for post-fringe
        # stop integration.
        post_fs_decorr = 1 - (1 - max_decorr) / (1 - pre_fs_decorr)
        int_time = np.sqrt(6 * post_fs_decorr / (np.pi ** 2 * max_rfac))

        # compute the number of samples that can be averaged using a power-of-two scheme
        num_two_foldings = int(
            np.floor(np.log2(int_time / corr_int_time.to(units.s).value))
        )
        return num_two_foldings
    else:
        # we're already above acceptable decorrelation value; cannot compress further
        return 0
