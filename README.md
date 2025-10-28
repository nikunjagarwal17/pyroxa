# PyroXa: Chemical Kinetics & Reactor Simulation

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![PyPI version](https://img.shields.io/pypi/v/pyroxa.svg)](https://pypi.org/project/pyroxa/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CI Tests](https://github.com/nikunjagarwal17/pyroxa/workflows/CI%20Tests/badge.svg)](https://github.com/nikunjagarwal17/pyroxa/actions)
[![Downloads](https://img.shields.io/pypi/dm/pyroxa.svg)](https://pypi.org/project/pyroxa/)

A pure Python library for chemical kinetics and reactor simulation. Built for students, researchers, and industry professionals. **Now available on PyPI!**

---

## 🚀 Quick Install

```bash
pip install pyroxa
```

Start using PyroXa in seconds - no compilation, no hassle!

---

## Features

- **Reaction Kinetics**: Elementary reactions, enzyme kinetics, surface reactions, temperature dependencies
- **Reactors**: Batch, CSTR, PFR, Packed Bed, Fluidized Bed, reactor networks
- **Thermodynamics**: Heat capacity, enthalpy, equilibrium, equations of state, phase equilibrium
- **Transport**: Diffusion, heat/mass transfer, dimensionless numbers
- **Analysis**: Sensitivity analysis, optimization, statistical methods, numerical solvers
- **Specialized**: Crystallization, catalyst performance, RTD, scale-up, process safety

132+ functions covering all aspects of chemical reaction engineering.

---

## Installation

### From PyPI (Recommended)

PyroXa is now available on PyPI! Install with a single command:

```bash
pip install pyroxa
```

That's it! You can now import and use PyroXa in your Python projects.

### Verify Installation

```bash
python -c "import pyroxa; print(f'PyroXa v{pyroxa.__version__} installed successfully!')"
```

### From Source (For Development)

If you want to contribute or modify the source code:

```bash
git clone https://github.com/nikunjagarwal17/pyroxa.git
cd pyroxa
pip install -e .
```

**Requirements**: Python 3.8+, NumPy, SciPy, Matplotlib (optional), PyYAML

See [INSTALLATION_GUIDE.md](./INSTALLATION_GUIDE.md) for detailed instructions and troubleshooting.

---

## Quick Start

### Installation
```bash
pip install pyroxa
```

### Basic Usage

```python
import pyroxa

# Rate constant using Arrhenius equation
k = pyroxa.arrhenius_rate(A=1e10, Ea=50000, T=298.15)
print(f"Rate constant: {k:.2e} 1/s")

# Batch reactor - calculate time for conversion
time = pyroxa.batch_reactor_time(initial_conc=1.0, final_conc=0.2, rate_constant=0.15, order=1)
print(f"Time required: {time:.2f} time units")

# CSTR - calculate volume for target conversion
volume = pyroxa.cstr_volume(flow_rate=10.0, rate_constant=0.5, conversion=0.8, order=1)
print(f"CSTR volume: {volume:.2f} volume units")

# PFR - calculate volume for target conversion
pfr_vol = pyroxa.pfr_volume(F_A0=10.0, X=0.8, C_A0=1.0, k=0.5, order=1)
print(f"PFR volume: {pfr_vol:.2f} volume units")

# Heat capacity using NASA polynomials
cp = pyroxa.heat_capacity_nasa(T=500.0, coeffs=[3.298677, 0.0014082404, -3.963222e-06, 5.641515e-09, -2.444854e-12])
print(f"Heat capacity: {cp:.2f} J/mol/K")
```

More examples in [`examples/`](./examples/) folder and on [PyPI](https://pypi.org/project/pyroxa/).

---

## Documentation

- **[Installation Guide](./INSTALLATION_GUIDE.md)** - Setup and troubleshooting
- **[Function Documentation](./DOCUMENTATION.md)** - Complete API reference
- **[Examples](./examples/)** - Industrial, pharmaceutical, and transport examples

---

## Testing

PyroXa includes a comprehensive test suite covering all 132+ functions and classes:

```bash
# Quick validation (fast)
python tests/quick_test.py

# Run all tests
python -m pytest tests/ -v

# Run specific category
python -m pytest tests/test_basic_kinetics.py -v
```

**Test Coverage:**
- 24 test files covering 26 documentation categories
- 200+ individual test methods
- 100% category coverage
- Core classes, exceptions, and utilities tested
- Automatic CI testing on Windows, macOS, Linux

See [tests/README.md](./tests/README.md) for complete test documentation.

---

## Project Structure

```
pyroxa/              # Main library (132 functions)
examples/            # Example scripts
tests/               # Test suite
docs/                # Documentation
```

---

## Use Cases

- **Students**: Homework, projects, learning reaction engineering concepts
- **Researchers**: Mechanism validation, parameter estimation, model development
- **Industry**: Reactor design, optimization, scale-up, safety analysis

**Get started in seconds:**
```bash
pip install pyroxa
python -c "import pyroxa; print(pyroxa.arrhenius_rate(A=1e10, Ea=50000, T=298.15))"
```

---

## Contributing

Contributions welcome! Fork the repo, make changes, add tests, and submit a PR.

---

## Authors

**Nikunj Agarwal** ([@nikunjagarwal17](https://github.com/nikunjagarwal17)) - Lead Developer  
**Contact**: [nikunjagarwal1704@gmail.com](mailto:nikunjagarwal1704@gmail.com)  
**PyPI**: [pypi.org/project/pyroxa](https://pypi.org/project/pyroxa/)  
**Issues**: [GitHub Issues](https://github.com/nikunjagarwal17/pyroxa/issues)

---

## Links

- **PyPI Package**: https://pypi.org/project/pyroxa/
- **GitHub Repository**: https://github.com/nikunjagarwal17/pyroxa
- **Documentation**: [DOCUMENTATION.md](./DOCUMENTATION.md)
- **Installation Guide**: [INSTALLATION_GUIDE.md](./INSTALLATION_GUIDE.md)

---

## License

MIT License - see LICENSE file for details.

---

**PyroXa v1.0.0** | Pure Python | Chemical Engineering Made Simple | `pip install pyroxa`
