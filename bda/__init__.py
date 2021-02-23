# -*- coding: utf-8 -*-
# Copyright (c) 2018 Paul La Plante
# Licensed under the 2-clause BSD License
"""Package for applying baseline-dependent averaging to radio astronomy datasets."""

__all__ = ["bda_tools", "decorr_calc", "apply_bda"]

from . import decorr_calc
from . import bda_tools
from .bda_tools import apply_bda
