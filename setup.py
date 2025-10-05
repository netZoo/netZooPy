from setuptools import setup, find_packages

setup(name='netZooPy',
    version='0.10.9',
    description='Python implementation of netZoo',
    url='https://github.com/netZoo/netZooPy',
    author='netZoo team',
    author_email='vfanfani@hsph.harvard.edu',
    license='GPL-3',
    long_description="""# netZooPy\n netZooPy is a python package to reconstruct, analyse, and plot biological networks.\n netZooPy currently integrates
    (gpu)PANDA, (gpu)LIONESS, (gpu)PUMA, SAMBAR, CONDOR, OTTER, DRAGON, COBRA, and BONOBO.\n""",
    long_description_content_type='text/markdown',
    packages=find_packages(),
    install_requires=[
    'h5py',    
    'pandas',
    'numpy>=1.19.2,<2.0',
    'networkx>=2.6.3',
    'matplotlib>=3.3.4',
    'scipy>=1.5.3',
    'igraph',
    'joblib>=1.1.0',
    'statsmodels>=0.12.2',
    'click',
    'tables'
    ],
    zip_safe=False,
    # add cli interface
    entry_points={
        'console_scripts': [
            'netzoopy=netZooPy.cli:cli'
        ],
    },   )
