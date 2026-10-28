"""Bind the primary water-monolayer parameter to the existing source inventory.

This changes initial parent coverage before depletion. It neither rescales
already solved fragments/spectra nor supplies chemistry, electron histories,
physical monolayer definitions or an experimental best-fit value.
"""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from .coverage_inventory import (
    InventoryPoolAllocation, MonolayerDefinition, ParentSurfaceInventory,
    build_coverage_source_sweep, prepare_coverage_inventory,
)

PARAMETER_NAME = 'initial_h2o_monolayers'


def bind_initial_h2o_monolayers(surface, initial_h2o_monolayers, *, monolayer,
        allocations, state_species, patch_multiplier=1., provenance,
        status='unvalidated', subdivisions_per_interval=1):
    """Prepare scalar or increasing nonnegative monolayer candidates.

    Every allocation must map to H2O. Fractions split one H2O parent inventory
    among water pools/internal states and are checked by the existing owner;
    they are not fragment yields. Actual patch areas enter exactly once.
    """
    if not isinstance(monolayer, MonolayerDefinition):
        raise ValueError('a separately defined MonolayerDefinition is required')
    allocations = tuple(allocations)
    if not allocations or any(not isinstance(a, InventoryPoolAllocation) for a in allocations):
        raise ValueError('supply explicit H2O pool allocations')
    if not isinstance(state_species, Mapping):
        raise ValueError('source states require an explicit state-to-species mapping')
    ids = {a.reservoir.state_id for a in allocations}
    if set(state_species) != ids or any(state_species[key] != 'H2O' for key in ids):
        raise ValueError('common-water source allocations must map exactly to H2O; fragments come from chemistry')
    parent = ParentSurfaceInventory(parent_id='H2O',monolayer=monolayer,
        allocations=allocations,patch_multiplier=patch_multiplier,
        profile_provenance=provenance,status=status,provenance=provenance)
    values = np.atleast_1d(initial_h2o_monolayers)
    return prepare_coverage_inventory(surface,parent,values,
        subdivisions_per_interval=subdivisions_per_interval)


def build_common_water_source_sweep(surface, initial_h2o_monolayers, *,
        parameter_provenance, **binding_inputs):
    """Recompute the actual depletion history at every monolayer candidate.

    The returned existing CoverageSourceSweep is the source input to complete
    transport/CR/optical cases. It is not a universal amplitude-scaling model.
    """
    plan=bind_initial_h2o_monolayers(surface,initial_h2o_monolayers,**binding_inputs)
    return build_coverage_source_sweep(plan.surface,plan.parent,plan.coverage_ml,
        parameter_provenance=parameter_provenance,
        subdivisions_per_interval=plan.subdivisions_per_interval)
