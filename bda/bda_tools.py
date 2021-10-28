# -*- coding: utf-8 -*-
# Copyright (c) 2019 Paul La Plante
# Licensed under the 2-clause BSD License
"""Main tools for applying BDA."""

import numpy as np
from astropy import units
from astropy.time import Time
import astropy.constants as const
from astropy.coordinates import Angle, EarthLocation, SkyCoord
from astropy.units import UnitConversionError

from pyuvdata import UVData
import pyuvdata.utils as uvutils

from . import decorr_calc as dc


def apply_bda(
    uv, max_decorr, pre_fs_int_time, corr_fov_angle, max_time, corr_int_time=None
):
    """Apply baseline dependent averaging to a UVData object.

    For each baseline in the UVData object, the expected decorrelation from
    averaging in time is computed. Baselines are averaged together in powers-
    of-two until the specified level of decorrelation is reached (rounded down).

    Parameters
    ----------
    uv : UVData object
        The UVData object to apply BDA to. No changes are made to this object,
        and instead a copy is returned.
    max_decorr : float
        The maximum decorrelation fraction desired in the output object. Must
        be between 0 and 1.
    pre_fs_int_time : astropy Quantity
        The pre-finge-stopping integration time inside of the correlator. The
        quantity should be compatible with units of time.
    corr_fov_angle : astropy Angle
        The opening angle at which the maximum decorrelation is to be
        calculated. Because a priori it is not known in which direction the
        decorrelation will be largest, the expected decorrelation is computed in
        all 4 cardinal directions at `corr_fov_angle` degrees off of zenith,
        and the largest one is used. This is a "worst case scenario"
        decorrelation.
    max_time : astropy Quantity
        The maximum amount of time that spectra from different times should be
        combined for. The ultimate integration time for a given baseline will be
        for max_time or the integration time that is smaller than the specified
        decorrelation level, whichever is smaller. The quantity should be
        compatible with units of time.
    corr_int_time : astropy Quantity, optional
        The output time of the correlator. If not specified, the smallest
        integration_time in the UVData object is used. If specified, the
        quantity should be compatible with units of time.

    Returns
    -------
    uv2 : UVData object
        The UVData object with BDA applied.

    Raises
    ------
    ValueError
        This is raised if the input parameters are not the appropriate type or
        in the appropriate range. It is also raised if the input UVData object
        is not in drift mode (the BDA code does rephasing within an averaged
        set of baselines).
    AssertionError
        This is raised if the baselines of the UVData object are not time-
        ordered.
    """
    if not isinstance(uv, UVData):
        raise ValueError(
            "apply_bda must be passed a UVData object as its first argument"
        )
    if not isinstance(corr_fov_angle, Angle):
        raise ValueError(
            "corr_fov_angle must be an Angle object from astropy.coordinates"
        )
    if not isinstance(pre_fs_int_time, units.Quantity):
        raise ValueError("pre_fs_int_time must be an astropy.units.Quantity")
    try:
        pre_fs_int_time.to(units.s)
    except UnitConversionError:
        raise ValueError("pre_fs_int_time must be a Quantity with units of time")
    if (
        corr_fov_angle.to(units.deg).value < 0
        or corr_fov_angle.to(units.deg).value > 90
    ):
        raise ValueError("corr_fov_angle must be between 0 and 90 degrees")
    if max_decorr < 0 or max_decorr > 1:
        raise ValueError("max_decorr must be between 0 and 1")
    if not isinstance(max_time, units.Quantity):
        raise ValueError("max_time must be an astropy.units.Quantity")
    try:
        max_time.to(units.s)
    except UnitConversionError:
        raise ValueError("max_time must be a Quantity with units of time")
    if corr_int_time is None:
        # assume the correlator integration time is the smallest int_time of the
        # UVData object
        corr_int_time = np.unique(uv.integration_time)[0] * units.s
    else:
        if not isinstance(corr_int_time, units.Quantity):
            raise ValueError("corr_int_time must be an astropy.units.Quantity")
        try:
            corr_int_time.to(units.s)
        except UnitConversionError:
            raise ValueError("corr_int_time must be a Quantity with units of time")
    if uv.phase_type != "drift":
        raise ValueError("UVData object must be in drift mode to apply BDA")

    # get relevant bits of metadata
    freq = np.amax(uv.freq_array[0, :]) * units.Hz
    chan_width = uv.channel_width * units.Hz
    antpos_enu, ants = uv.get_ENU_antpos()
    lat, lon, alt = uv.telescope_location_lat_lon_alt
    antpos_ecef = uvutils.ECEF_from_ENU(antpos_enu, lat, lon, alt)
    telescope_location = EarthLocation.from_geocentric(
        uv.telescope_location[0],
        uv.telescope_location[1],
        uv.telescope_location[2],
        unit="m",
    )

    # make a new UVData object to put BDA baselines in
    uv2 = UVData()

    # copy over metadata
    uv2.Nbls = uv.Nbls
    uv2.Nfreqs = uv.Nfreqs
    uv2.Npols = uv.Npols
    uv2.vis_units = uv.vis_units
    uv2.Nspws = uv.Nspws
    uv2.spw_array = uv.spw_array
    uv2.freq_array = uv.freq_array
    uv2.polarization_array = uv.polarization_array
    uv2.channel_width = uv.channel_width
    uv2.object_name = uv.object_name
    uv2.telescope_name = uv.telescope_name
    uv2.instrument = uv.instrument
    uv2.telescope_location = uv.telescope_location
    history = uv.history + " Baseline dependent averaging applied."
    uv2.history = history
    uv2.Nants_data = uv.Nants_data
    uv2.Nants_telescope = uv.Nants_telescope
    uv2.antenna_names = uv.antenna_names
    uv2.antenna_numbers = uv.antenna_numbers
    uv2.x_orientation = uv.x_orientation
    uv2.extra_keywords = uv.extra_keywords
    uv2.antenna_positions = uv.antenna_positions
    uv2.antenna_diameters = uv.antenna_diameters
    uv2.gst0 = uv.gst0
    uv2.rdate = uv.rdate
    uv2.earth_omega = uv.earth_omega
    uv2.dut1 = uv.dut1
    uv2.timesys = uv.timesys
    uv2.uvplane_reference_time = uv.uvplane_reference_time

    # initialize place-keeping variables and Nblt-sized metadata
    start_index = 0
    uv2.Nblts = 0
    uv2.uvw_array = np.zeros_like(uv.uvw_array)
    uv2.time_array = np.zeros_like(uv.time_array)
    uv2.lst_array = np.zeros_like(uv.lst_array)
    uv2.ant_1_array = np.zeros_like(uv.ant_1_array)
    uv2.ant_2_array = np.zeros_like(uv.ant_2_array)
    uv2.baseline_array = np.zeros_like(uv.baseline_array)
    uv2.integration_time = np.zeros_like(uv.integration_time)
    uv2.data_array = np.zeros_like(uv.data_array)
    uv2.flag_array = np.zeros_like(uv.flag_array)
    uv2.nsample_array = np.zeros_like(uv.nsample_array)

    # iterate over baselines
    for key in uv.get_antpairs():
        print("averaging baseline ", key)
        ind1, ind2, indp = uv._key2inds(key)
        if len(ind2) != 0:
            raise AssertionError(
                "ind2 from _key2inds() is not 0--exiting. This should not happen, "
                "please contact the package maintainers."
            )
        data = uv._smart_slicing(
            uv.data_array, ind1, ind2, indp, squeeze="none", force_copy=True
        )
        flags = uv._smart_slicing(
            uv.flag_array, ind1, ind2, indp, squeeze="none", force_copy=True
        )
        nsamples = uv._smart_slicing(
            uv.nsample_array, ind1, ind2, indp, squeeze="none", force_copy=True
        )

        # get lx and ly for baseline
        ant1 = np.where(ants == key[0])[0][0]
        ant2 = np.where(ants == key[1])[0][0]
        x1, y1, z1 = antpos_ecef[ant1, :]
        x2, y2, z2 = antpos_ecef[ant2, :]
        lx = np.abs(x2 - x1) * units.m
        ly = np.abs(y2 - y1) * units.m

        # figure out how many time samples we can combine together
        if key[0] == key[1]:
            # autocorrelation--don't average
            n_two_foldings = 0
        else:
            n_two_foldings = dc.bda_compression_factor(
                max_decorr,
                freq,
                lx,
                ly,
                corr_fov_angle,
                chan_width,
                pre_fs_int_time,
                corr_int_time,
            )
        # convert from max_time to max_samples
        max_samples = (max_time / corr_int_time).to(units.dimensionless_unscaled)
        max_two_foldings = int(np.floor(np.log2(max_samples)))
        n_two_foldings = min(n_two_foldings, max_two_foldings)
        n_int = 2 ** (n_two_foldings)
        print("averaging {:d} time samples...".format(n_int))

        # figure out how many output samples we're going to have
        n_in = len(ind1)
        n_out = n_in // n_int + min(1, n_in % n_int)

        # get relevant metdata
        uvw_array = uv.uvw_array[ind1, :]
        times = uv.time_array[ind1]
        if not np.all(times == np.sort(times)):
            raise AssertionError(
                "times of uvdata object are not monotonically increasing; "
                "throwing our hands up"
            )
        lsts = uv.lst_array[ind1]
        int_time = uv.integration_time[ind1]

        # do the averaging
        input_shape = data.shape
        assert input_shape == (n_in, 1, uv.Nfreqs, uv.Npols)
        output_shape = (n_out, 1, uv.Nfreqs, uv.Npols)
        data_out = np.empty(output_shape, dtype=np.complex128)
        flags_out = np.empty(output_shape, dtype=np.bool_)
        nsamples_out = np.empty(output_shape, dtype=np.float32)
        uvws_out = np.empty((n_out, 3), dtype=np.float64)
        times_out = np.empty((n_out,), dtype=np.float64)
        lst_out = np.empty((n_out,), dtype=np.float64)
        int_time_out = np.empty((n_out,), dtype=np.float64)

        if n_out == n_in:
            # we don't need to average
            current_index = start_index + n_out
            uv2.data_array[start_index:current_index, :, :, :] = data
            uv2.flag_array[start_index:current_index, :, :, :] = flags
            uv2.nsample_array[start_index:current_index, :, :, :] = nsamples
            uv2.uvw_array[start_index:current_index, :] = uvw_array
            uv2.time_array[start_index:current_index] = times
            uv2.lst_array[start_index:current_index] = lsts
            uv2.integration_time[start_index:current_index] = int_time
            uv2.ant_1_array[start_index:current_index] = key[0]
            uv2.ant_2_array[start_index:current_index] = key[1]
            uv2.baseline_array[start_index:current_index] = uvutils.antnums_to_baseline(
                ant1, ant2, None
            )
            start_index = current_index

        else:
            # rats, we actually have to do work...

            # phase up the data along each chunk of times
            for i in range(n_out):
                # compute zenith of the desired output time
                i1 = i * n_int
                i2 = min((i + 1) * n_int, n_in)
                assert i2 - i1 > 0
                t0 = Time((times[i1] + times[i2 - 1]) / 2, scale="utc", format="jd")
                zenith_coord = SkyCoord(
                    alt=Angle(90 * units.deg),
                    az=Angle(0 * units.deg),
                    obstime=t0,
                    frame="altaz",
                    location=telescope_location,
                )
                obs_zenith_coord = zenith_coord.transform_to("icrs")
                zenith_ra = obs_zenith_coord.ra
                zenith_dec = obs_zenith_coord.dec

                # get data, flags, and nsamples of slices
                data_chunk = data[i1:i2, :, :, :]
                flags_chunk = flags[i1:i2, :, :, :]
                nsamples_chunk = nsamples[i1:i2, :, :, :]

                # actually phase now
                # compute new uvw coordinates
                icrs_coord = SkyCoord(
                    ra=zenith_ra, dec=zenith_dec, unit="radian", frame="icrs"
                )
                uvws = np.float64(uvw_array[i1:i2, :])
                itrs_telescope_location = SkyCoord(
                    x=uv.telescope_location[0] * units.m,
                    y=uv.telescope_location[1] * units.m,
                    z=uv.telescope_location[2] * units.m,
                    representation_type="cartesian",
                    frame="itrs",
                    obstime=t0,
                )
                itrs_lat_lon_alt = uv.telescope_location_lat_lon_alt

                frame_telescope_location = itrs_telescope_location.transform_to("icrs")

                frame_telescope_location.representation_type = "cartesian"

                uvw_ecef = uvutils.ECEF_from_ENU(uvws, *itrs_lat_lon_alt)

                itrs_uvw_coord = SkyCoord(
                    x=uvw_ecef[:, 0] * units.m,
                    y=uvw_ecef[:, 1] * units.m,
                    z=uvw_ecef[:, 2] * units.m,
                    representation_type="cartesian",
                    frame="itrs",
                    obstime=t0,
                )
                frame_uvw_coord = itrs_uvw_coord.transform_to("icrs")

                frame_rel_uvw = (
                    frame_uvw_coord.cartesian.get_xyz().value.T
                    - frame_telescope_location.cartesian.get_xyz().value
                )

                new_uvws = uvutils.phase_uvw(
                    icrs_coord.ra.rad, icrs_coord.dec.rad, frame_rel_uvw
                )

                # average these uvws together to get the "average" position in
                # the uv-plane
                avg_uvws = np.average(new_uvws, axis=0)

                # calculate and apply phasor
                w_lambda = (
                    new_uvws[:, 2].reshape((i2 - i1), 1)
                    / const.c.to("m/s").value
                    * uv.freq_array.reshape(1, uv.Nfreqs)
                )
                phs = np.exp(-1j * 2 * np.pi * w_lambda[:, None, :, None])
                data_chunk *= phs

                # sum data, propagate flag array, and adjusting nsample accordingly
                data_slice = np.sum(data_chunk, axis=0)
                flag_slice = np.sum(flags_chunk, axis=0)
                nsamples_slice = np.sum(nsamples_chunk, axis=0) / (i2 - i1)
                data_out[i, :, :, :] = data_slice
                flags_out[i, :, :, :] = flag_slice
                nsamples_out[i, :, :, :] = nsamples_slice

                # update metadata
                uvws_out[i, :] = avg_uvws
                times_out[i] = (times[i1] + times[i2 - 1]) / 2
                lst_out[i] = (lsts[i1] + lsts[i2 - 1]) / 2
                int_time_out[i] = np.average(int_time[i1:i2]) * (i2 - i1)

            # update data and metadata when we're done with this baseline
            current_index = start_index + n_out
            uv2.data_array[start_index:current_index, :, :, :] = data_out
            uv2.flag_array[start_index:current_index, :, :, :] = flags_out
            uv2.nsample_array[start_index:current_index, :, :, :] = nsamples_out
            uv2.uvw_array[start_index:current_index, :] = uvws_out
            uv2.time_array[start_index:current_index] = times_out
            uv2.lst_array[start_index:current_index] = lst_out
            uv2.integration_time[start_index:current_index] = int_time_out
            uv2.ant_1_array[start_index:current_index] = key[0]
            uv2.ant_2_array[start_index:current_index] = key[1]
            uv2.baseline_array[start_index:current_index] = uvutils.antnums_to_baseline(
                ant1, ant2, None
            )
            start_index = current_index

    # clean up -- shorten all arrays to actually be size nblts
    nblts = start_index
    uv2.Nblts = nblts
    uv2.data_array = uv2.data_array[:nblts, :, :, :]
    uv2.flag_array = uv2.flag_array[:nblts, :, :, :]
    uv2.nsample_array = uv2.nsample_array[:nblts, :, :, :]
    uv2.uvw_array = uv2.uvw_array[:nblts, :]
    uv2.time_array = uv2.time_array[:nblts]
    uv2.lst_array = uv2.lst_array[:nblts]
    uv2.integration_time = uv2.integration_time[:nblts]
    uv2.ant_1_array = uv2.ant_1_array[:nblts]
    uv2.ant_2_array = uv2.ant_2_array[:nblts]
    uv2.baseline_array = uv2.baseline_array[:nblts]
    uv2.Ntimes = len(np.unique(uv2.time_array))

    # set phasing info
    uv2.phase_type = "phased"
    uv2.phase_center_ra = zenith_ra.rad
    uv2.phase_center_dec = zenith_dec.rad
    uv2.phase_center_frame = 2000.0

    # fix up to correct old phasing method
    uv2.phase(zenith_ra.rad, zenith_dec.rad, epoch="J2000", fix_old_proj=True)

    # run a check
    uv2.check()

    return uv2
