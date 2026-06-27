# Boundary-case report

Version 2 classifies boundary cases instead of forcing every edge or ranking into a binary decision.

## Edge statuses

- `retained`
- `boundary`
- `rejected`

## Node and ranking statuses

- `resolved`
- `boundary`
- `unresolved`
- `bridge_like`

## Boundary reasons

The implementation writes explicit `boundary_reason` values such as:

- `low_rank_margin`
- `threshold_ci_overlap`
- `unstable_representation_support`
- `mixed_component_membership`
- `weak_weighted_consensus`

## Scientific status

This is a development classification layer requiring user inspection.

