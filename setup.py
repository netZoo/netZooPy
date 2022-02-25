from setuptools import setup, find_packages

setup(name='netZooPy',
    version='0.9.2',
    description='Python implementation of netZoo.',
    url='https://github.com/netZoo/netZooPy',
    author='netZoo team',
    author_email='vfanfani@hsph.harvard.edu',
    license='GPL-3',
    packages=find_packages(),
    install_requires=['pandas',
    'numpy',
    'networkx',
    'matplotlib',
    'scipy',
    'python-igraph',
    'joblib',
    'statsmodels',
    'click'
    ],
    zip_safe=False,
    # add cli interface
    entry_points={
        'console_scripts': [
            'netzoopy=netZooPy.cli:cli'
        ],
    },   )
