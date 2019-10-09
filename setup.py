from setuptools import setup, find_packages

setup(name='netZooPy',
    version='0.1.1',
    description='Python implementation of PUMA, PANDA, LIONESS.',
    url='https://github.com/netZoo/netZooPy',
    author='Alessandro Marin, Cho-Yi Chen, David van IJzendoorn',
    author_email='twangxx@hsph.harvard.edu',
    license='GPL-3',
    packages=['netZooPy'],
    install_requires=['pandas',
    'numpy',
    'networkx',
    'matplotlib',
    'scipy'
    ],
    zip_safe=False)
