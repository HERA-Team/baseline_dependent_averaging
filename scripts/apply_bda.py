#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2018 Paul La Plante
# Licensed under the 2-clause BSD License

from __future__ import print_function, division, absolute_import

import argparse
import os
import sys
from pyuvdata import UVData
from astropy import units
from astropy.coordinates import Angle

from bda import bda_tools


# setup argparse
a = argparse.ArgumentParser(
    description="A command-line script for performing baseline-dependent averaging"
)
a.add_argument("--file_in", type=str, help="input data file")
a.add_argument("--file_out", type=str, help="output data file")
a.add_argument(
    "--max_decorr",
    type=float,
    default=0.1,
    help="maximum amount of decorrelation allowed; default is 10%",
)
a.add_argument(
    "--pre_fs_int_time",
    type=float,
    default=0.1,
    help="time in seconds of phase stopping in correlator",
)
a.add_argument(
    "--corr_FoV_angle",
    type=float,
    default=20.0,
    help="FoV angle in degrees at which to compute max_decorr; default is 20",
)
a.add_argument(
    "--max_time",
    type=float,
    default=30.0,
    help="maximum amount of time to average; default is 30 seconds",
)
a.add_argument(
    "--corr_int_time",
    type=float,
    default=None,
    required=False,
    help="total integration time of correlator; defaults to smallest integration_time in file",
)
a.add_argument(
    "--overwrite",
    default=False,
    action="store_true",
    help="overwrite output file if it already exists",
)
a.add_argument(
    "--filetype",
    default="uvh5",
    type=str,
    help="output filetype, should be 'uvh5', 'uvfits', or 'miriad'",
)

# get args
args = a.parse_args()

if os.path.exists(args.file_out) and args.overwrite is False:
    print("{} exists. Use --overwrite to overwrite the file.".format(args.file_out))
    sys.exit(0)

# check that output filetype is valid
if args.filetype not in ("uvh5", "uvfits", "miriad"):
    print("filetype must be one of uvh5, uvfits, or miriad")
    sys.exit(0)

# read in file
uv = UVData()
uv.read(args.file_in)

# apply BDA
pre_fs_int_time = args.pre_fs_int_time * units.s
corr_FoV_angle = Angle(args.corr_FoV_angle, units.deg)
max_time = args.max_time * units.s
uv2 = bda_tools.apply_bda(
    uv, args.max_decorr, pre_fs_int_time, corr_FoV_angle, max_time, args.corr_int_time
)

# write out file
if args.filetype == "uvh5":
    uv2.write_uvh5(args.file_out, clobber=True)
if args.filetype == "uvfits":
    uv2.write_uvfits(args.file_out, spoof_nonessential=True, force_phase=True)
if args.filetype == "miriad":
    uv2.write_miriad(args.file_out, clobber=True)
