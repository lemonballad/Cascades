"""Core computation modules for spectroscopy simulations."""

from cascades.core.basis import rcbasis, basis_tc
from cascades.core.franck_condon import fcfac2_tc, fcinfo_tc
from cascades.core.response import cascade_2drr_res, LaserParameters, MaterialParameters
from cascades.core.fsrs import cascade_fsrs_res, FSRSLaserParameters, FSRSMaterialParameters
from cascades.core.offres import cascade_2drr_offres, OffResMaterialParameters
from cascades.core.optimized import fcfac2_tc_fast, fcinfo_tc_fast, check_numba_available

__all__ = [
    "rcbasis", "basis_tc",
    "fcfac2_tc", "fcinfo_tc",
    "fcfac2_tc_fast", "fcinfo_tc_fast", "check_numba_available",
    "cascade_2drr_res", "LaserParameters", "MaterialParameters",
    "cascade_fsrs_res", "FSRSLaserParameters", "FSRSMaterialParameters",
    "cascade_2drr_offres", "OffResMaterialParameters",
]
