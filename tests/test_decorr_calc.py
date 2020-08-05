# -*- coding: utf-8 -*-
# Copyright (c) 2020 The HERA Collaboration
# Licensed under the 2-clause BSD License

import pytest
import numpy as np
from astropy import units, constants
from astropy.coordinates import Angle

from bda import decorr_calc as dc


def test_dudt():
    # define quantities
    lx = 14.0 * units.m
    ly = 14.0 * units.m
    hour_angle = Angle(0.0, unit="rad")
    earth_omega = 2 * np.pi * units.rad / units.sday
    wavelength = constants.c / (250 * units.MHz)

    # compute reference value
    dudt_check = (
        lx.value * earth_omega.to("rad/s").value / wavelength.to("m").value
    )

    # run function and check
    dudt = dc._dudt(lx, ly, hour_angle, earth_omega, wavelength)
    assert isinstance(dudt, units.Quantity)
    assert np.isclose(dudt.to("rad/s").value, dudt_check)

    return


def test_dvdt():
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
