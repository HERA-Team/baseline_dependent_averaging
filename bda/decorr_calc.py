# -*- coding: utf-8 -*-
# Copyright (c) 2018 Paul La Plante
# Licensed under the 2-clause BSD License

from __future__ import print_function, division, absolute_import

import numpy as np
from astropy import constants as const
from astropy.coordinates import Angle
from astropy import units

# define some HERA-specific constants
hera_latitude = Angle("-30:43:17.5", units.deg)


def _dudt(lx, ly, hour_angle, earth_omega, wavelength):
    """
    Dosctring goes here
    """
    hour_angle = max(
        min(hour_angle, Angle(90.0, units.degree)), Angle(-90.0, units.degree)
    )
    return (
        (lx * np.cos(hour_angle) - ly * np.sin(hour_angle)) * earth_omega / wavelength
    )


def _dvdt(lx, ly, hour_angle, dec, earth_omega, wavelength):
    """
    Docstring goes here
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
    """
    Docstring goes here
    """
    wavelength = const.c / frequency.to(1 / units.s)
    earth_rot_speed = (Angle(360, units.deg) / units.sday).to(units.arcminute / units.s)
    max_resolution = Angle(np.arcsin(wavelength / baseline), units.rad)
    return float(pre_fs_int_time * earth_rot_speed / max_resolution.to(units.arcminute))


def decorr_chan_width(chan_width, baseline, corr_FoV):
    """
    Docstring goes here
    """
    return float(
        chan_width.to(1 / units.s) * baseline * np.sin(corr_FoV.to(units.rad)) / const.c
    )


def decorr_post_fs_int_time(
    lx, ly, post_fs_int_time, corr_FoV, frequency, telescope_latitude=hera_latitude
):
    """
    Docstring goes here
    """
    wavelength = const.c / frequency.to(1 / units.s)
    earth_rot_speed = (Angle(360, units.deg) / units.sday).to(units.arcminute / units.s)

    # case 1: +l
    du = _dudt(lx, ly, corr_FoV, earth_rot_speed, wavelength)
    l = np.cos(90.0 * units.deg + corr_FoV)
    rfac = (du * l) ** 2

    # case 2: -l
    du = _dudt(lx, ly, -corr_FoV, earth_rot_speed, wavelength)
    l = np.cos(90.0 * units.deg - corr_FoV)
    rfac = max(rfac, (du * l) ** 2)

    # case 3: +m
    dv = _dvdt(lx, ly, 0.0, telescope_latitude + corr_FoV, earth_rot_speed, wavelength)
    m = np.cos(90.0 * units.deg + corr_FoV)
    rfac = max(rfac, (dv * m) ** 2)

    # case 4: -m
    dv = _dvdt(lx, ly, 0.0, telescope_latitude - corr_FoV, earth_rot_speed, wavelength)
    m = np.cos(90.0 * units.deg - corr_FoV)
    rfac = max(rfac, (dv * m) ** 2)

    # make sure we have the right units
    rfac = rfac.to(units.rad ** 2 / units.s ** 2)

    # add other factors; return time and max rfac value
    decorr_frac = np.pi ** 2 * (post_fs_int_time.to(units.s).value) ** 2 / 6.0 * rfac.value
    return decorr_frac, rfac.value


def bda_compression_factor(
    max_decorr=0.1,
    frequency=(250.0 * 1e6 * units.Hz),
    lx=(14.6 * units.m),
    ly=(14.6 * units.m),
    corr_FoV_angle=Angle(20.0, units.degree),
    chan_width=(30.517 * units.kHz),
    pre_fs_int_time=(0.1 * units.s),
    corr_int_time=(10 * units.s),
):

    # calculate the pre-BDA decorrelation given the correlator settings
    baseline = np.sqrt(lx ** 2 + ly ** 2)
    decorr_cw = decorr_chan_width(chan_width, baseline, corr_FoV_angle)

    decorr_pre_int = decorr_pre_fs_int_time(frequency, baseline, pre_fs_int_time)

    decorr_post_int, max_rfac = decorr_post_fs_int_time(
        lx, ly, corr_int_time, corr_FoV_angle, frequency
    )

    # calculate pre- and post-fs decorrelations
    pre_fs_decorr = 1 - (1 - decorr_cw) * (1 - decorr_pre_int)
    total_decorr = 1 - (1 - pre_fs_decorr) * (1 - decorr_post_int)

    if total_decorr < max_decorr:
        # if total decorrelation is less than max allowed, we can average
        # figure out the maximum amount of decorrelation allowed for post-fringe stop integration
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
