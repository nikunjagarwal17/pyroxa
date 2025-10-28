"""
Quick Test - PyroXa v1.0.0
Fast validation of core functionality
"""

import sys
import os
import traceback
import warnings

# Suppress NumPy warnings for Python 3.13
warnings.filterwarnings('ignore', category=RuntimeWarning, module='numpy')
warnings.filterwarnings('ignore', message='.*MINGW.*')

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

def test_import():
    """Test PyroXa import"""
    try:
        import pyroxa
        version = getattr(pyroxa, '__version__', '1.0.0')
        print(f"✓ PyroXa v{version} imported successfully")
        return True
    except Exception as e:
        print(f"✗ Import failed: {e}")
        traceback.print_exc()
        return False

def test_function_count():
    """Test that functions are available"""
    try:
        import pyroxa
        # Count all non-private attributes
        all_items = dir(pyroxa)
        functions = [item for item in all_items if not item.startswith('_')]
        count = len(functions)
        print(f"✓ {count} functions/classes available")
        return True
    except Exception as e:
        print(f"✗ Function count test failed: {e}")
        traceback.print_exc()
        return False

def test_arrhenius_rate():
    """Test Arrhenius rate calculation"""
    try:
        import pyroxa
        k = pyroxa.arrhenius_rate(A=1e10, Ea=50000, T=298.15)
        if k > 0:
            print(f"✓ Arrhenius rate: {k:.4e} 1/s")
            return True
        else:
            print(f"✗ Invalid arrhenius rate: {k}")
            return False
    except Exception as e:
        print(f"✗ Arrhenius rate test failed: {e}")
        traceback.print_exc()
        return False

def test_batch_reactor():
    """Test batch reactor time calculation"""
    try:
        import pyroxa
        time = pyroxa.batch_reactor_time(initial_conc=1.0, final_conc=0.2, rate_constant=0.1, order=1)
        if time > 0:
            print(f"✓ Batch reactor time: {time:.2f} s for 80% conversion")
            return True
        else:
            print(f"✗ Invalid time: {time}")
            return False
    except Exception as e:
        print(f"✗ Batch reactor test failed: {e}")
        traceback.print_exc()
        return False

def test_cstr_volume():
    """Test CSTR volume calculation"""
    try:
        import pyroxa
        volume = pyroxa.cstr_volume(flow_rate=0.1, rate_constant=0.1, conversion=0.8, order=1)
        if volume > 0:
            print(f"✓ CSTR volume: {volume:.3f} L for 80% conversion")
            return True
        else:
            print(f"✗ Invalid volume: {volume}")
            return False
    except Exception as e:
        print(f"✗ CSTR volume test failed: {e}")
        traceback.print_exc()
        return False

def test_heat_capacity():
    """Test thermodynamic calculation"""
    try:
        import pyroxa
        coeffs = [3.298677, 0.0014082404, -3.963222e-06, 5.641515e-09, -2.444854e-12, 0.0, 0.0]
        cp = pyroxa.heat_capacity_nasa(T=500.0, coeffs=coeffs)
        if cp > 0:
            print(f"✓ Heat capacity: {cp:.2f} J/(mol·K)")
            return True
        else:
            print(f"✗ Invalid heat capacity: {cp}")
            return False
    except Exception as e:
        print(f"✗ Heat capacity test failed: {e}")
        traceback.print_exc()
        return False

def test_reynolds_number():
    """Test Reynolds number calculation"""
    try:
        import pyroxa
        Re = pyroxa.reynolds_number(density=1000, velocity=1.0, length=0.1, viscosity=0.001)
        if Re > 0:
            print(f"✓ Reynolds number: {Re:.0f}")
            return True
        else:
            print(f"✗ Invalid Reynolds number: {Re}")
            return False
    except Exception as e:
        print(f"✗ Reynolds number test failed: {e}")
        traceback.print_exc()
        return False

def main():
    """Run quick tests"""
    print("=" * 60)
    print("PyroXa Quick Test Suite")
    print("=" * 60)
    
    tests = [
        ("Import Test", test_import),
        ("Function Count", test_function_count),
        ("Arrhenius Rate", test_arrhenius_rate),
        ("Batch Reactor Time", test_batch_reactor),
        ("CSTR Volume", test_cstr_volume),
        ("Heat Capacity", test_heat_capacity),
        ("Reynolds Number", test_reynolds_number),
    ]
    
    results = []
    for name, test_func in tests:
        print(f"\n{name}:")
        try:
            results.append(test_func())
        except Exception as e:
            print(f"✗ Test crashed: {e}")
            traceback.print_exc()
            results.append(False)
    
    print("\n" + "=" * 60)
    passed = sum(results)
    total = len(results)
    print(f"Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("✓ All tests passed!")
        print("=" * 60)
        return 0
    else:
        print(f"✗ {total - passed} test(s) failed")
        print("=" * 60)
        return 1

if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        print(f"\n✗ FATAL ERROR: {e}")
        traceback.print_exc()
        sys.exit(1)
