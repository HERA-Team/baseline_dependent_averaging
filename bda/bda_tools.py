# -*- coding: utf-8 -*-
# Copyright (c) 2018 Paul La Plante
# Licensed under the 2-clause BSD License

from __future__ import print_function, division, absolute_import

import numpy as np
from astropy import units
from astropy.time import Time
import astropy.constants as const
from astropy.coordinates import Angle, EarthLocation, SkyCoord
from pyuvdata import UVData
import pyuvdata.utils as uvutils

from . import decorr_calc as dc


def apply_bda(uv, max_decorr, pre_fs_int_time, corr_FoV_angle, max_samples, corr_int_time=None):
    """
    Apply baseline dependent averaging to a UVData object.
    """
    if not isinstance(uv, UVData):
        raise ValueError("apply_bda must be passed a UVData object as its first argument")
    if not isinstance(corr_FoV_angle, Angle):
        raise ValueError("corr_FoV_angle must be an Angle object from astropy.coordinates")
    if not isinstance(pre_fs_int_time, units.Quantity):
        raise ValueError("pre_fs_int_time must be an astropy.units.Quantity")
    try:
        pre_fs_int_time.to(units.s)
    except UnitConversionError:
        raise ValueError("pre_fs_int_time must be a Quantity with units of time")
    if corr_FoV_angle.to(units.deg).value < 0 or corr_FoV_angle.to(units.deg).value > 90:
        raise ValueError("corr_FoV_angle must be between 0 and 90 degrees")
    if max_decorr < 0 or max_decorr > 1:
        raise ValueError("max_decorr must be between 0 and 1")
    if corr_int_time is None:
        # assume the correlator integration time is the smallest int_time of the UVData object
        corr_int_time = np.unique(uv.integration_time)[0] * units.s
    else:
        if not isinstance(corr_int_time, units.Quantity):
            raise ValueError("corr_int_time must be an astropy.units.Quantity")
        try:
            corr_int_time.to(units.s)
        except UnitConversionError:
            raise ValueError("corr_int_time must be a Quantity with units of time")

    # get relevant bits of metadata
    freq = np.amax(uv.freq_array[0, :]) * units.Hz
    chan_width = uv.channel_width * units.Hz
    antpos_enu, ants = uv.get_ENU_antpos()
    telescope_location = EarthLocation.from_geocentric(uv.telescope_location[0],
                                                       uv.telescope_location[1],
                                                       uv.telescope_location[2],
                                                       unit='m')

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
    history = uv.history + ' Baseline dependent averaging applied.'
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
    uv2.uvw_array = uv.uvw_array
    uv2.time_array = uv.time_array
    uv2.lst_array = uv.lst_array
    uv2.ant_1_array = uv.ant_1_array
    uv2.ant_2_array = uv.ant_2_array
    uv2.baseline_array = uv.baseline_array
    uv2.integration_time = uv.integration_time
    uv2.data_array = uv.data_array
    uv2.flag_array = uv.flag_array
    uv2.nsample_array = uv.nsample_array

    # iterate over baselines
    for key in uv.get_antpairs():
        print("averaging baseline ", key)
        ind1, ind2, indp = uv._key2inds(key)
        assert len(ind2) == 0
        data = uv._smart_slicing(uv.data_array, ind1, ind2, indp, squeeze='none', force_copy=True)
        flags = uv._smart_slicing(uv.flag_array, ind1, ind2, indp, squeeze='none', force_copy=True)
        nsamples = uv._smart_slicing(uv.nsample_array, ind1, ind2, indp, squeeze='none', force_copy=True)

        # get lx and ly for baseline
        ant1 = np.where(ants==key[0])[0][0]
        ant2 = np.where(ants==key[1])[0][0]
        x1, y1, z1 = antpos_enu[ant1, :]
        x2, y2, z2 = antpos_enu[ant2, :]
        lx = np.abs(x2 - x1) * units.m
        ly = np.abs(y2 - y1) * units.m

        # figure out how many time samples we can combine together
        if key[0] == key[1]:
            # autocorrelation--don't average
            n_int = 1
        else:
            n_int = dc.bda_compression_factor(max_decorr, freq, lx, ly, corr_FoV_angle, chan_width,
                                              pre_fs_int_time, corr_int_time)
        n_int = min(n_int, max_samples)
        print("averaging {:d} time samples...".format(n_int))

        # figure out how many output samples we're going to have
        n_in = len(ind1)
        n_out = n_in // n_int + min(1, n_in % n_int)

        # get relevant metdata
        uvw_array = uv.uvw_array[ind1, :]
        times = uv.time_array[ind1]
        if not np.all(times == np.sort(times)):
            raise AssertionError("times of uvdata object are not monotonically increasing; "
                                 "throwing our hands up")
        lsts = uv.lst_array[ind1]
        int_time = uv.integration_time[ind1]
        baselines = uv.baseline_array[ind1]

        # do the averaging
        input_shape = data.shape
        assert input_shape == (n_in, 1, uv.Nfreqs, uv.Npols)
        output_shape = (n_out, 1, uv.Nfreqs, uv.Npols)
        data_out = np.empty(output_shape, dtype=np.complex128)
        flags_out = np.empty(output_shape, dtype=np.bool)
        nsamples_out = np.empty(output_shape, dtype=np.float32)
        uvws_out = np.empty((n_out, 3), dtype=np.float64)
        times_out = np.empty((n_out,), dtype=np.float64)
        lst_out = np.empty((n_out,), dtype=np.float64)
        int_time_out = np.empty((n_out,), dtype=np.float64)

        # phase up the data along each chunk of times
        for i in range(n_out):
            # compute zenith of the desired output time
            i1 = i * n_int
            i2 = min((i + 1) * n_int, n_in)
            t0 = Time((times[i1] + times[i2 - 1]) / 2, scale='utc', format='jd')
            zenith_coord = SkyCoord(alt=Angle(90 * units.deg), az=Angle(0 * units.deg),
                                    obstime=t0, frame='altaz', location=telescope_location)
            obs_zenith_coord = zenith_coord.transform_to('icrs')
            zenith_ra = obs_zenith_coord.ra
            zenith_dec = obs_zenith_coord.dec

            # get data, flags, and nsamples of slices
            data_chunk = data[i1:i2, :, :, :]
            flags_chunk = flags[i1:i2, :, :, :]
            nsamples_chunk = nsamples[i1:i2, :, :, :]

            # actually phase now
            # compute new uvw coordinates
            icrs_coord = SkyCoord(ra=zenith_ra, dec=zenith_dec, unit='radian', frame='icrs')
            uvws = np.float64(uvw_array[i1:i2, :])
            itrs_telescope_location = SkyCoord(x=uv.telescope_location[0] * units.m,
                                               y=uv.telescope_location[1] * units.m,
                                               z=uv.telescope_location[2] * units.m,
                                               representation='cartesian',
                                               frame='itrs', obstime=t0)
            itrs_lat_lon_alt = uv.telescope_location_lat_lon_alt

            frame_telescope_location = itrs_telescope_location.transform_to('icrs')

            frame_telescope_location.representation = 'cartesian'

            uvw_ecef = uvutils.ECEF_from_ENU(uvws, *itrs_lat_lon_alt)

            itrs_uvw_coord = SkyCoord(x=uvw_ecef[:, 0] * units.m,
                                      y=uvw_ecef[:, 1] * units.m,
                                      z=uvw_ecef[:, 2] * units.m,
                                      representation='cartesian',
                                      frame='itrs', obstime=t0)
            frame_uvw_coord = itrs_uvw_coord.transform_to('icrs')

            frame_rel_uvw = (frame_uvw_coord.cartesian.get_xyz().value.T
                             - frame_telescope_location.cartesian.get_xyz().value)

            new_uvws = uvutils.phase_uvw(icrs_coord.ra.rad,
                                         icrs_coord.dec.rad,
                                         frame_rel_uvw)

            # average these uvws together to get the "average" position in the uv-plane
            avg_uvws = np.average(new_uvws, axis=0)

            # calculate and apply phasor
            w_lambda = (new_uvws[:, 2].reshape((i2 - i1), 1)
                        / const.c.to('m/s').value * uv.freq_array.reshape(1, uv.Nfreqs))
            phs = np.exp(-1j * 2 * np.pi * w_lambda[:, None, :, None])
            data_chunk *= phs

            # sum data, propagate flag array, and adjusting nsample accordingly
            data_slice = np.average(data_chunk, axis=0)
            flag_slice = np.sum(flags_chunk, axis=0)
            nsamples_slice = np.average(nsamples_chunk, axis=0)
            data_out[i, :, :, :] = data_slice
            flags_out[i, :, :, :] = flag_slice
            nsamples_out[i, :, :, :] = nsamples_slice

            # update metadata
            uvws_out[i, :] = avg_uvws
            times_out[i] = (times[i1] + times[i2 - 1]) / 2
            lst_out[i] = (lsts[i1] + lsts[i2 - 1]) / 2
            int_time_out[i] = np.average(int_time[i1:i2]) * (i2 - i1 - 1)

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
        uv2.baseline_array[start_index:current_index] = uvutils.antnums_to_baseline(ant1, ant2,
                                                                                    None)
        start_index = current_index

    # clean up -- shorten all arrays to actually be size Nblts
    Nblts = start_index
    uv2.Nblts = Nblts
    uv2.data_array = uv2.data_array[:Nblts, :, :, :]
    uv2.flag_array = uv2.flag_array[:Nblts, :, :, :]
    uv2.nsample_array = uv2.nsample_array[:Nblts, :, :, :]
    uv2.uvw_array = uv2.uvw_array[:Nblts, :]
    uv2.time_array = uv2.time_array[:Nblts]
    uv2.lst_array = uv2.lst_array[:Nblts]
    uv2.integration_time = uv2.integration_time[:Nblts]
    uv2.ant_1_array = uv2.ant_1_array[:Nblts]
    uv2.ant_2_array = uv2.ant_2_array[:Nblts]
    uv2.baseline_array = uv2.baseline_array[:Nblts]
    uv2.Ntimes = len(np.unique(uv2.time_array))

    # run a check
    uv2.check()

    return uv2
