# -*- coding: utf-8 -*-
# Copyright (c) 2021 The HERA Collaboration
# Licensed under the 2-clause BSD License
"""Tests for decorr_calc.py module."""

import numpy as np
from astropy import units, constants
from astropy.coordinates import Angle

from bda import decorr_calc as dc


def test_dudt():
    """Test _dudt function."""
    # define quantities
    lx = 14.0 * units.m
    ly = 14.0 * units.m
    hour_angle = Angle(0.0, unit="rad")
    earth_omega = 2 * np.pi * units.rad / units.sday
    wavelength = constants.c / (250 * units.MHz)

    # compute reference value
    dudt_check = lx.value * earth_omega.to("rad/s").value / wavelength.to("m").value

    # run function and check
    dudt = dc._dudt(lx, ly, hour_angle, earth_omega, wavelength)
    assert isinstance(dudt, units.Quantity)
    assert np.isclose(dudt.to("rad/s").value, dudt_check)

    return


def test_dvdt():
    """Test _dvdt function."""
    # define quantities
    lx = 14.0 * units.m
    ly = 14.0 * units.m
    hour_angle = Angle(0.0, unit="rad")
    dec = Angle(-30.0, unit="deg")
    earth_omega = 2 * np.pi * units.rad / units.sday
    wavelength = constants.c / (250 * units.MHz)

    # compute reference value
    dvdt_check = (
        -ly.value * 0.5 * earth_omega.to("rad/s").value / wavelength.to("m").value
    )

    # run function and check
    dvdt = dc._dvdt(lx, ly, hour_angle, dec, earth_omega, wavelength)
    assert isinstance(dvdt, units.Quantity)
    assert np.isclose(dvdt.to("rad/s").value, dvdt_check)

    return


def test_decorr_pre_fs_int_time():
    """Test decorr_pre_fs_int_time function."""
    # define quantities
    freq = 250 * units.MHz
    bl = 14.0 * units.m
    pre_fs_int_time = 1 * units.s
    decorr = dc.decorr_pre_fs_int_time(freq, bl, pre_fs_int_time)

    # compute comparison value
    earth_omega = (2 * np.pi * units.rad / units.sday).to(units.arcminute / units.s)
    wavelength = constants.c / freq.to(1 / units.s)
    max_res = Angle(np.arcsin(wavelength / bl), units.rad)
    decorr_ref = float(pre_fs_int_time * earth_omega / max_res.to(units.arcminute))
    assert np.isclose(decorr, decorr_ref)

    return


def test_decorr_chan_width():
    """Test decorr_chan_width function."""
    # define quantities
    corr_fov = Angle(20.0, units.deg)
    bl = 14.0 * units.m
    chan_width = 50 * units.kHz
    decorr = dc.decorr_chan_width(chan_width, bl, corr_fov)

    # compute comparison value
    decorr_ref = (
        chan_width.to(1 / units.s) * bl * np.sin(corr_fov.to(units.rad)) / constants.c
    )
    assert np.isclose(decorr, float(decorr_ref))

    return


def test_decorr_post_fs_int_time():
    """Test decorr_post_fs_int_time function."""
    # define quantities
    lx = 14.0 * units.m
    ly = 14.0 * units.m
    post_fs_int_time = 10 * units.s
    corr_fov = Angle(20.0, units.deg)
    frequency = 250 * units.MHz
    decorr_frac, rfac = dc.decorr_post_fs_int_time(
        lx, ly, post_fs_int_time, corr_fov, frequency
    )

    # compute comparison value
    wavelength = constants.c / frequency.to(1 / units.s)
    earth_rot_speed = (Angle(360, units.deg) / units.sday).to(units.arcminute / units.s)
    du = dc._dudt(lx, ly, -corr_fov, earth_rot_speed, wavelength)
    lval = np.cos(90.0 * units.deg - corr_fov)
    rfac_ref = ((du * lval) ** 2).to(units.rad ** 2 / units.s ** 2).value
    assert np.isclose(rfac, rfac_ref)

    decorr_frac_ref = (
        np.pi ** 2 * (post_fs_int_time.to(units.s).value) ** 2 / 6.0 * rfac_ref
    )
    assert np.isclose(decorr_frac_ref, decorr_frac)

    return


def test_bda_compression_factor():
    """Test bda_compression_factor function."""
    # define quantities
    max_decorr = 0.1
    frequency = 250 * units.MHz
    lx = 14.0 * units.m
    ly = 14.0 * units.m
    corr_fov = Angle(20.0, units.deg)
    chan_width = 50 * units.kHz
    pre_fs_int_time = 0.1 * units.s
    corr_int_time = 10 * units.s
    num_two_foldings = dc.bda_compression_factor(
        max_decorr,
        frequency,
        lx,
        ly,
        corr_fov,
        chan_width,
        pre_fs_int_time,
        corr_int_time,
    )
    assert num_two_foldings == 6

    # now test a long baseline that can't be averaged more
    lx = 1 * units.km
    num_two_foldings = dc.bda_compression_factor(
        max_decorr,
        frequency,
        lx,
        ly,
        corr_fov,
        chan_width,
        pre_fs_int_time,
        corr_int_time,
    )
    assert num_two_foldings == 0

    return
