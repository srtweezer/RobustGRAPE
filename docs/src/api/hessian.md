# Hessian Analysis

```@meta
CurrentModule = RobustGRAPE
```

Near an optimum the fidelity landscape is low rank: for typical Rydberg Hamiltonians only
a handful of control directions change the fidelity at all. The functions on this page
expose those directions, which is what makes a closed-loop calibration tractable — instead
of searching ``\mathbb{R}^{2L}`` for a length-``L`` pulse with amplitude and phase control,
one searches the few leading eigenvectors.

Writing ``\bm{u}`` for a small deviation from parameters that give ``F = 1``, the fidelity
Hessian ``\mathcal{H} = -\partial^2 F/\partial\bm{u}\,\partial\bm{u}`` gives

```math
F(\bm{u}) \approx 1 - \tfrac{1}{2}\, \bm{u}^T \mathcal{H}\, \bm{u}.
```

For a gate whose point is *robustness* the relevant quadratic form is larger. With an error
``\varepsilon \sim \mathcal{N}(0, \sigma_\varepsilon^2)`` the average fidelity is
``\langle F\rangle = F - \sigma_\varepsilon^2 S_\varepsilon``, so a deviation costs

```math
1 - \langle F(\bm{u})\rangle = \tfrac{1}{2}\,\bm{u}^T \mathcal{H}\, \bm{u}
                             + \tfrac{1}{2}\,\sigma_\varepsilon^2\, \bm{u}^T \mathcal{Q}\, \bm{u},
```

where ``\mathcal{Q} = \partial^2 S_\varepsilon/\partial\bm{u}\,\partial\bm{u}`` is the
**robustness Hessian**. The principal directions of ``\mathcal{H}`` are the ones a
calibration should search; those of ``\mathcal{Q}`` are the ones it must avoid, or it
spends the robustness the gate time was paid for. Because the two do not share
eigenvectors, the total form ``\mathcal{H} + \sigma_\varepsilon^2 \mathcal{Q}`` has to be
diagonalized at the specific ``\sigma_\varepsilon`` of interest.

## Fidelity Hessian

```@docs
calculate_fidelity_hessian
calculate_principal_parameters
```

## Robustness Hessian and principal waveforms

```@docs
calculate_robustness_hessian
principal_waveforms
hessian_health
```

## Nuisance parameters and reparameterization

A gate's free single-qubit phases, gauge angles and target-family parameters are
*nuisance* directions: they are re-optimized rather than calibrated, so they must be
eliminated from the Hessian rather than diagonalized with it. `schur_eliminate` does that
elimination; `optimize_additional_parameters` re-maximizes them first, and
`shift_error_source` evaluates the form at a displaced error so a robustness Hessian can
be built by differencing.

`apply_control_constraint` handles the case where the physical knobs are fewer than the
propagation's parameters — a composite gate that plays one waveform in each of several
segments propagates `n_segments × ntimes` steps but programs only `ntimes`.

```@docs
schur_eliminate
optimize_additional_parameters
shift_error_source
apply_control_constraint
```
