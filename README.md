# RobustGRAPE.jl

A Julia package for efficient quantum optimal control with error robustness analysis and optimization.

## Overview

RobustGRAPE.jl provides tools for designing high-fidelity quantum gates that are robust against noise and experimental errors. It is based on the GRadient Ascent Pulse Engineering (GRAPE) pulse algorithm: time-evolution is approximated by piecewise constant evolution operators. In addition to the bare GRAPE algorithm, this package offers tools to analyze the sensitivity of quantum gates to various errors, and can be used to design quantum gates with reduced sensitivity to certain errors.

This package was devised for quantum gates based on the Rydberg interaction between two atoms, but should be adaptable to many different kinds of quantum systems. It relies on analytical calculations to provide fast and accurate calculations of, e.g., the average gate fidelity and its sensitivity to error pulses, as well as the gradient of these values with respect to control parameters. On the other hand, the derivatives of each piecewise constant evolution operator are calculated using finite differences to ensure maximum adaptability to various systems. Therefore, although this package is meant to offer resonable performance, it is certainly possible to be much faster using code tailored to specific scenarios.

## References

This paper introduces the GRAPE algorithm: N. Khaneja, T. Reiss, C. Kehlet, T. Schulte-Herbrüggen, and S. J. Glaser, Optimal control of coupled spin dynamics: Design of NMR pulse sequences by gradient ascent algorithms, J. Magn. Reson. 172, 296 (2005).

If you use the fidelity response functions in your work, please cite: Richard Bing-Shiun Tsai, Xiangkai Sun, Adam L. Shaw, Ran Finkelstein, and Manuel Endres. “Benchmarking and Fidelity Response Theory of High-Fidelity Rydberg Entangling Gates.” PRX Quantum 6, no. 1 (2025): 010331. https://doi.org/10.1103/PRXQuantum.6.010331.

## Functionalities

### Fidelity calculations

The package provides a simple interface to compute:
- The fidelity of a parametrized quantum gate.
- The second-order sensitivity of the fidelity to a given noise operator.
- The gradient of these two values w.r.t control parameters.

Based on these features, the package provides a high-level interface to minimize a global cost that includes the infidelity and the sensitivity to error sources.

Additionally, the package offers error analysis tools:
- The average expectation value of a given operator (e.g., a projector onto a state that experiences some decay)
- The fidelity response function of the fidelity to a noise described by a classical power spectral density, such as a laser noise. We invite the reader to familiarize themselves with this tool in the above [reference](https://doi.org/10.1103/PRXQuantum.6.010331).

### Unitary calculations

At a lower level, the fidelity tools rely on the ability to efficiently compute the total unitary gate, the derivative of this unitary w.r.t. control parameters, the derivative of this unitary w.r.t. error sources, and the derivative of this unitary w.r.t. error sources and control parameters. For users that wish to work with a different quantity than the fidelity, they can directly access the corresponding methods.

### Regularization

Provides simple functionalities to regularize control parameters and their derivatives. Regularization can promote convergence to smooth high-fidelity pulses.

### Rydberg tools

Provides pre-defined Hamiltonian and parametrized CZ gates for a two-atom Rydberg system.

## Installation

This package is currently not available on the Julia Register, but can be installed directly by the user:

```julia
using Pkg
Pkg.add(url="https://github.com/srtweezer/RobustGRAPE")
```


## Quick start & Documentation

Please see the [online documentation](https://srtweezer.github.io/RobustGRAPE/).

For quick start, please consult the [Rydberg 2-atom time-optimal gate example](https://srtweezer.github.io/RobustGRAPE/dev/examples/), which showcases most features in the package.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

This package was developed with assistance from [Claude](https://claude.ai), particularly for documentation and code organization.