from setuptools import setup, find_packages

setup(name='netZooPy',
    version='0.5',
    description='Python implementation of netZoo.',
    url='https://github.com/netZoo/netZooPy',
    author='netZoo team',
    author_email='twangxx@hsph.harvard.edu',
    license='GPL-3',
    packages=['netZooPy'],
    install_requires=['pandas',
    'numpy',
    'networkx',
    'matplotlib',
    'scipy',
    'python-igraph'
    ],
    zip_safe=False)
