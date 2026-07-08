from setuptools import setup, find_packages

setup(
    name='icetemp',
    version='0.1.0',
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    install_requires=[
        'numpy',
        'pandas',
        'scipy',
        'matplotlib',
        'cmcrameri',
    ],
    extras_require={
        # Tier-3 Bayesian (Kennedy-O'Hagan) calibration, src/icetemp/calibration/.
        # Kept as an extra (not a hard install_requires) because surmise pins
        # numpy<2.2/scipy<1.15 -- narrower than the rest of icetemp needs.
        'calibration': [
            'scikit-learn>=1.2',
            'emcee>=3.1',
            'surmise>=0.4',
            'arviz>=0.23',
            'corner>=2.3',
        ],
    },
)
