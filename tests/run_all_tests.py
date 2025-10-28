"""
Run All PyroXa Tests
"""

import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def run_all_tests():
    """Run all test suites"""
    try:
        import pytest
    except ImportError:
        print("pytest not installed. Installing...")
        os.system("pip install pytest")
        import pytest
    
    # Run pytest with verbose output
    args = [
        '-v',  # Verbose
        '--tb=short',  # Short traceback
        '-x',  # Stop on first failure (optional, remove for full run)
        'tests/'
    ]
    
    return pytest.main(args)

if __name__ == "__main__":
    sys.exit(run_all_tests())
