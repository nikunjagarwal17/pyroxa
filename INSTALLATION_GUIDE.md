# PyroXa Installation Guide

Quick guide to get PyroXa installed on your system.

---

## Requirements

- **Python**: 3.8 or newer
- **Operating System**: Windows 10/11, macOS 10.14+, Linux (Ubuntu, Debian, CentOS)
- **RAM**: 2GB minimum (4GB+ recommended for large simulations)
- **Disk Space**: ~500 MB

Check your Python version:
```bash
python --version
```

---

## Quick Installation

```bash
# Clone repository
git clone https://github.com/nikunjagarwal17/pyroxa.git
cd pyroxa/project

# Install dependencies and PyroXa
pip install -r requirements.txt
pip install -e .

# Verify installation
python -c "import pyroxa; print(f'PyroXa v{pyroxa.__version__} installed successfully')"
```

---

## Virtual Environment (Recommended)

Using a virtual environment keeps PyroXa isolated from other Python projects.

**Using venv:**
```bash
# Create environment
python -m venv pyroxa_env

# Activate
# Windows:
pyroxa_env\Scripts\activate
# macOS/Linux:
source pyroxa_env/bin/activate

# Install PyroXa
pip install -r requirements.txt
pip install -e .

# Deactivate when done
deactivate
```

**Using conda:**
```bash
conda create -n pyroxa python=3.10
conda activate pyroxa
pip install -r requirements.txt
pip install -e .
```

---

## Step-by-Step Installation

### 1. Install Python
If Python is not installed or version < 3.8:
- **Windows**: Download from [python.org](https://www.python.org/downloads/)
- **Linux**: `sudo apt install python3 python3-pip`
- **macOS**: `brew install python3`

### 2. Clone Repository
```bash
git clone https://github.com/nikunjagarwal17/pyroxa.git
cd pyroxa/project
```

Or download ZIP from GitHub and extract.

### 3. Install Dependencies
```bash
pip install -r requirements.txt
```

### 4. Install PyroXa

**Recommended: Development/Editable Mode**
```bash
pip install -e .
```
This uses `setup.py` to install PyroXa in editable mode. Changes to source files take effect immediately without reinstalling.

**Alternative: Standard Installation**
```bash
pip install .
```
Installs PyroXa as a regular package using `setup.py`. Requires reinstallation after code changes.

**What does setup.py do?**
- Defines package metadata (name, version, dependencies)
- Specifies which files to include
- Handles installation process
- Ensures dependencies (NumPy, SciPy, etc.) are installed

### 5. Verify Installation
```bash
python -c "import pyroxa; print(pyroxa.__version__)"
```

---

## Testing

```bash
# Quick test
python tests/quick_test.py

# Full test suite
python -m pytest tests/ -v
```

---

## Troubleshooting

**"Module not found: pyroxa"**
```bash
cd pyroxa/project
pip install -e .
```

**Permission denied (macOS/Linux)**
```bash
pip install --user -r requirements.txt
pip install --user -e .
```

**NumPy/SciPy installation fails**
```bash
pip install --upgrade pip setuptools wheel
pip install -r requirements.txt
```

**DLL load failed (Windows)**
- Install [Microsoft Visual C++ Redistributable](https://aka.ms/vs/17/release/vc_redist.x64.exe)
- Restart computer
- Reinstall PyroXa

---

## Updating PyroXa

```bash
cd pyroxa/project
git pull origin main
pip install -e . --upgrade
```

---

## Uninstallation

```bash
pip uninstall pyroxa
```

---

## Getting Help

**Email**: [nikunjagarwal1704@gmail.com](mailto:nikunjagarwal1704@gmail.com)  
**Issues**: [GitHub Issues](https://github.com/nikunjagarwal17/pyroxa/issues)

Include:
- Operating system
- Python version
- Error message

---

## Next Steps

1. Try example: `python -c "import pyroxa; k = pyroxa.arrhenius_rate(1e10, 50000, 298.15); print(f'k = {k:.2e}')"`
2. Run examples: `python examples/industrial_showcase.py`
3. Read [README.md](./README.md) for overview
4. Check [FUNCTION_DOCUMENTATION.md](./FUNCTION_DOCUMENTATION.md) for API reference

---

**Nikunj Agarwal** | [nikunjagarwal1704@gmail.com](mailto:nikunjagarwal1704@gmail.com)
