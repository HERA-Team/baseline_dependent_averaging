# -*- coding: utf-8 -*-
# Copyright (c) 2021 The HERA Collaboration
# Licensed under the 2-clause BSD License
"""Tests for bda_tools.py module."""

import pytest
import numpy as np
from numpy.random import default_rng
from astropy import units
from astropy.coordinates import Angle
import pyuvdata
import pyuvdata.utils as uvutils
from pyuvdata import UVData

from bda import bda_tools


@pytest.fixture(scope="module")
def fake_data_generator():
    """Generate a fake dataset as a module-level fixture."""
    # generate a fake data file for testing BDA application
    uvd = UVData()
    t0 = 2459000.0
    dt = (2 * units.s).to(units.day)
    t1 = t0 + dt.value
    t2 = t0 + 2 * dt.value
    # define baseline-time spacing
    uvd.ant_1_array = np.asarray([0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=np.int_)
    uvd.ant_2_array = np.asarray([0, 1, 2, 0, 1, 2, 0, 1, 2], dtype=np.int_)
    uvd.time_array = np.asarray([t0, t0, t0, t1, t1, t1, t2, t2, t2], dtype=np.float64)
    uvd.Ntimes = 3
    uvd.Nbls = 3
    uvd.Nblts = uvd.Ntimes * uvd.Nbls

    # define frequency array
    nfreqs = 1024
    freq_array = np.linspace(50e6, 250e6, num=nfreqs)
    uvd.freq_array = freq_array[np.newaxis, :]
    uvd.spw_array = [0]
    uvd.Nfreqs = nfreqs
    uvd.Nspws = 1

    # define polarization array
    uvd.polarization_array = [-5]
    uvd.Npols = 1

    # make random data for data array
    data_shape = (uvd.Nblts, 1, uvd.Nfreqs, uvd.Npols)
    uvd.data_array = np.zeros(data_shape, dtype=np.complex128)
    rng = default_rng()
    uvd.data_array += rng.standard_normal(data_shape)
    uvd.data_array += 1j * rng.standard_normal(data_shape)
    uvd.flag_array = np.zeros_like(uvd.data_array, dtype=np.bool_)
    uvd.nsample_array = np.ones_like(uvd.data_array, dtype=np.float32)

    # set telescope and antenna positions
    hera_telescope = pyuvdata.telescopes.get_telescope("HERA")
    uvd.telescope_location_lat_lon_alt = hera_telescope.telescope_location_lat_lon_alt
    antpos0 = np.asarray([0, 0, 0], dtype=np.float64)
    antpos1 = np.asarray([14, 0, 0], dtype=np.float64)
    antpos2 = np.asarray([28, 0, 0], dtype=np.float64)
    antpos_enu = np.vstack((antpos0, antpos1, antpos2))
    antpos_xyz = uvutils.ECEF_from_ENU(antpos_enu, *uvd.telescope_location_lat_lon_alt)
    uvd.antenna_positions = antpos_xyz - uvd.telescope_location

    uvw_array = np.zeros((uvd.Nblts, 3), dtype=np.float64)
    uvw_array[::3, :] = antpos0 - antpos0
    uvw_array[1::3, :] = antpos1 - antpos0
    uvw_array[2::3, :] = antpos2 - antpos0
    uvd.uvw_array = uvw_array
    uvd.antenna_numbers = np.asarray([0, 1, 2], dtype=np.int_)
    uvd.antenna_names = np.asarray(["H0", "H1", "H2"], dtype=np.str_)
    uvd.Nants_data = 3
    uvd.Nants_telescope = 3

    # set other metadata
    uvd.vis_units = "uncalib"
    uvd.channel_width = 5e4  # 50 kHz
    uvd.phase_type = "drift"
    uvd.baseline_array = uvd.antnums_to_baseline(uvd.ant_1_array, uvd.ant_2_array)
    uvd.history = "BDA test file"
    uvd.instrument = "HERA"
    uvd.telescope_name = "HERA"
    uvd.object_name = "zenith"
    uvd.integration_time = 2 * np.ones_like(uvd.baseline_array, dtype=np.float64)
    uvd.set_lsts_from_time_array()

    # run a check
    uvd.check()

    yield uvd

    # clean up when done
    del uvd

    return


@pytest.fixture(scope="function")
def fake_data(fake_data_generator):
    """Make a per-function copy of main fake dataset."""
    uvd = fake_data_generator.copy()
    yield uvd
    del uvd
    return


def test_apply_bda(fake_data):
    """Generate test data and make sure that BDA is applied correctly."""
    uvd = fake_data

    # define parameters
    max_decorr = 0.1
    pre_fs_int_time = 0.1 * units.s
    corr_fov_angle = Angle(20.0, units.deg)
    max_time = 30 * units.s
    uvd2 = bda_tools.apply_bda(
        uvd,
        max_decorr,
        pre_fs_int_time,
        corr_fov_angle,
        max_time,
    )

    # make sure that things are correct
    assert uvd2.Ntimes == 3
    assert uvd2.Nbls == 3
    assert uvd2.Nblts == 5  # baselines (0, 1) and (0, 2) were all averaged together
    bl1_inds, _, _ = uvd._key2inds((0, 1))
    bl2_inds, _, _ = uvd2._key2inds((0, 1))
    target = uvd2.data_array[bl2_inds]
    # multiply in factor of 3 to account for 3 time samples being averaged
    avg_data = 3.0 * np.average(uvd.data_array[bl1_inds], axis=0)[np.newaxis, :, :, :]
    # there are small differences to account for phasing considerations, but the
    # results are still "close" for the times and baseline lengths chosen
    assert np.allclose(target, avg_data, atol=1e-3)

    return


def test_apply_bda_non_uvd_error():
    """Test error for not using a UVData object."""
    # define parameters
    uvd = "foo"
    max_decorr = 0.1
    pre_fs_int_time = 0.1 * units.s
    corr_fov_angle = Angle(20.0, units.degree)
    max_time = 30 * units.s
    corr_int_time = 2 * units.s

    # test using something besides a UVData object
    with pytest.raises(ValueError) as cm:
        bda_tools.apply_bda(
            uvd,
            max_decorr,
            pre_fs_int_time,
            corr_fov_angle,
            max_time,
            corr_int_time,
        )
    assert str(cm.value).startswith(
        "apply_bda must be passed a UVData object as its first argument"
    )
    return


def test_apply_bda_non_angle_error():
    """Test error for not using an Angle for corr_fov_angle."""
    # define parameters
    uvd = UVData()
    max_decorr = 0.1
    pre_fs_int_time = 0.1 * units.s
    corr_fov_angle = "foo"
    max_time = 30 * units.s
    corr_int_time = 2 * units.s

    # test using something besides an angle for corr_fov_angle
    with pytest.raises(ValueError) as cm:
        bda_tools.apply_bda(
            uvd,
            max_decorr,
            pre_fs_int_time,
            corr_fov_angle,
            max_time,
            corr_int_time,
        )
    assert str(cm.value).startswith(
        "corr_fov_angle must be an Angle object from astropy.coordinates"
    )
    return


def test_apply_bda_pre_fs_int_time_float_error():
    """Test error for not using a Quantity for pre_fs_int_time."""
    # define parameters
    uvd = UVData()
    max_decorr = 0.1
    pre_fs_int_time = 0.1
    corr_fov_angle = Angle(20.0, units.degree)
    max_time = 30 * units.s
    corr_int_time = 2 * units.s

    # pass in pre_fs_int_time as a float instead of a Quantity
    with pytest.raises(ValueError) as cm:
        bda_tools.apply_bda(
            uvd,
            max_decorr,
            pre_fs_int_time,
            corr_fov_angle,
            max_time,
            corr_int_time,
        )
    assert str(cm.value).startswith("pre_fs_int_time must be an astropy.units.Quantity")
    return


def test_apply_bda_pre_fs_int_time_bad_quantity_error():
    """Test error for using an incompatible Quantity for pre_fs_int_time."""
    # define parameters
    uvd = UVData()
    max_decorr = 0.1
    pre_fs_int_time = 0.1 * units.m
    corr_fov_angle = Angle(20.0, units.degree)
    max_time = 30 * units.s
    corr_int_time = 2 * units.s

    # use pre_fs_int_time as a quantity with bogus units
    with pytest.raises(ValueError) as cm:
        bda_tools.apply_bda(
            uvd,
            max_decorr,
            pre_fs_int_time,
            corr_fov_angle,
            max_time,
            corr_int_time,
        )
    assert str(cm.value).startswith(
        "pre_fs_int_time must be a Quantity with units of time"
    )
    return


def test_apply_bda_bad_corr_fov_angle_error():
    """Test error for using an invalid corr_fov_angle."""
    # define parameters
    uvd = UVData()
    max_decorr = 0.1
    pre_fs_int_time = 0.1 * units.s
    corr_fov_angle = Angle(180.0, units.degree)
    max_time = 30 * units.s
    corr_int_time = 2 * units.s

    # use a bad corr_fov_angle value
    with pytest.raises(ValueError) as cm:
        bda_tools.apply_bda(
            uvd,
            max_decorr,
            pre_fs_int_time,
            corr_fov_angle,
            max_time,
            corr_int_time,
        )
    assert str(cm.value).startswith("corr_fov_angle must be between 0 and 90 degrees")
    return


def test_apply_bda_max_decorr_error():
    """Test error for supplying an invalid max_decorr value."""
    # define parameters
    uvd = UVData()
    max_decorr = -0.1
    pre_fs_int_time = 0.1 * units.s
    corr_fov_angle = Angle(20.0, units.degree)
    max_time = 30 * units.s
    corr_int_time = 2 * units.s

    # use a bad max_decorr value
    with pytest.raises(ValueError) as cm:
        bda_tools.apply_bda(
            uvd,
            max_decorr,
            pre_fs_int_time,
            corr_fov_angle,
            max_time,
            corr_int_time,
        )
    assert str(cm.value).startswith("max_decorr must be between 0 and 1")
    return


def test_apply_bda_max_time_float_error():
    """Test error for using max_time as a float."""
    # define parameters
    uvd = UVData()
    max_decorr = 0.1
    pre_fs_int_time = 0.1 * units.s
    corr_fov_angle = Angle(20.0, units.degree)
    max_time = 30.0
    corr_int_time = 2 * units.s

    # use a plain float for max_time instead of a Quantity
    with pytest.raises(ValueError) as cm:
        bda_tools.apply_bda(
            uvd,
            max_decorr,
            pre_fs_int_time,
            corr_fov_angle,
            max_time,
            corr_int_time,
        )
    assert str(cm.value).startswith("max_time must be an astropy.units.Quantity")
    return


def test_apply_bda_max_time_bad_quantity_error():
    """Test error for using an incompatible quantity for max_time."""
    # define parameters
    uvd = UVData()
    max_decorr = 0.1
    pre_fs_int_time = 0.1 * units.s
    corr_fov_angle = Angle(20.0, units.degree)
    max_time = 30 * units.m
    corr_int_time = 2 * units.s

    # use a bad unit for max_time
    with pytest.raises(ValueError) as cm:
        bda_tools.apply_bda(
            uvd,
            max_decorr,
            pre_fs_int_time,
            corr_fov_angle,
            max_time,
            corr_int_time,
        )
    assert str(cm.value).startswith("max_time must be a Quantity with units of time")
    return


def test_apply_bda_corr_int_time_float_error():
    """Test error for using corr_int_time as a float."""
    # define parameters
    uvd = UVData()
    max_decorr = 0.1
    pre_fs_int_time = 0.1 * units.s
    corr_fov_angle = Angle(20.0, units.degree)
    max_time = 30 * units.s
    corr_int_time = 2.0

    # use a plain float for corr_int_time
    with pytest.raises(ValueError) as cm:
        bda_tools.apply_bda(
            uvd,
            max_decorr,
            pre_fs_int_time,
            corr_fov_angle,
            max_time,
            corr_int_time,
        )
    assert str(cm.value).startswith("corr_int_time must be an astropy.units.Quantity")
    return


def test_apply_bda_corr_int_time_bad_quantity_error():
    """Test error for using an incompatible quantity for corr_int_time."""
    # define parameters
    uvd = UVData()
    max_decorr = 0.1
    pre_fs_int_time = 0.1 * units.s
    corr_fov_angle = Angle(20.0, units.degree)
    max_time = 30 * units.s
    corr_int_time = 2.0 * units.m

    # use a bad quantity for corr_int_time
    with pytest.raises(ValueError) as cm:
        bda_tools.apply_bda(
            uvd,
            max_decorr,
            pre_fs_int_time,
            corr_fov_angle,
            max_time,
            corr_int_time,
        )
    assert str(cm.value).startswith(
        "corr_int_time must be a Quantity with units of time"
    )
    return


def test_apply_bda_phased_error(fake_data):
    """Test error when applying bda to phased data."""
    # convert input data to phased
    uvd = fake_data
    uvd.phase_to_time(uvd.time_array[0])
    max_decorr = 0.1
    pre_fs_int_time = 0.1 * units.s
    corr_fov_angle = Angle(20.0, units.degree)
    max_time = 30 * units.s
    corr_int_time = 2.0 * units.s

    with pytest.raises(ValueError) as cm:
        bda_tools.apply_bda(
            uvd,
            max_decorr,
            pre_fs_int_time,
            corr_fov_angle,
            max_time,
            corr_int_time,
        )
    assert str(cm.value).startswith("UVData object must be in drift mode to apply BDA")
    return


def test_apply_bda_ind2_key_error(fake_data):
    """Test error when the length of ind2 is non-zero."""
    # mess with fake data
    uvd = fake_data
    uvd.ant_1_array[1] = 1
    uvd.ant_2_array[1] = 0
    max_decorr = 0.1
    pre_fs_int_time = 0.1 * units.s
    corr_fov_angle = Angle(20.0, units.degree)
    max_time = 30 * units.s
    corr_int_time = 2.0 * units.s

    with pytest.raises(AssertionError) as cm:
        bda_tools.apply_bda(
            uvd,
            max_decorr,
            pre_fs_int_time,
            corr_fov_angle,
            max_time,
            corr_int_time,
        )
    assert str(cm.value).startswith("ind2 from _key2inds() is not 0--exiting")
    return


def test_apply_bda_non_increasing_error(fake_data):
    """Test error when times in object are not monotonically increasing."""
    # mess with fake data
    uvd = fake_data
    uvd.time_array = uvd.time_array[::-1]
    max_decorr = 0.1
    pre_fs_int_time = 0.1 * units.s
    corr_fov_angle = Angle(20.0, units.degree)
    max_time = 30 * units.s

    with pytest.raises(AssertionError) as cm:
        bda_tools.apply_bda(
            uvd,
            max_decorr,
            pre_fs_int_time,
            corr_fov_angle,
            max_time,
        )
    assert str(cm.value).startswith(
        "times of uvdata object are not monotonically increasing"
    )
    return
