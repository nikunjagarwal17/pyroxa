from setuptools import setup, find_packages
from setuptools.command.install import install

class PostInstallCommand(install):
    """Post-installation message"""
    def run(self):
        install.run(self)
        print("\nInstallation complete!")

# Pure Python PyroXa - Chemical Kinetics Library
# No compilation required - easy installation on any platform!

install_requires = [
    'numpy>=1.19.0',
    'scipy>=1.7.0',  # Added for advanced numerical methods
    'matplotlib>=3.3.0',  # Added for plotting capabilities
]

# Read README for long description
try:
    with open("README.md", "r", encoding="utf-8") as fh:
        long_description = fh.read()
except FileNotFoundError:
    long_description = "PyroXa: Pure Python chemical kinetics and reactor simulation library"

setup(
    name='pyroxa',
    version='1.0.0',
    packages=find_packages(),
    install_requires=install_requires,
    description='Pure Python chemical kinetics and reactor simulation library',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='Nikunj Agarwal',
    author_email='nikunjagarwal1704@gmail.com',
    url='https://github.com/nikunjagarwal17/pyroxa',
    python_requires='>=3.8',
    zip_safe=True,
    include_package_data=True,
    package_data={
        'pyroxa': ['*.py'],
    },
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Education',
    ],
    keywords='chemical kinetics, reactor simulation, chemical engineering, thermodynamics, pure python',
    cmdclass={
        'install': PostInstallCommand,
    },
)
