# Representation weighting report

The Version 2 workflow estimates non-negative representation weights constrained to sum to one.

## Current implementation

The current implementation uses a deterministic simplex-normalised pilot objective:

```text
w_r >= 0
sum(w_r) = 1
S_weighted(i,j) = sum_r w_r * S_r(i,j)
```

The objective rewards neighbourhood stability and cancer-type agreement, and penalises representation redundancy, batch association, and instability.

## Scientific status

This is a pilot approximation until a reviewed convex optimisation formulation and solver are approved.

