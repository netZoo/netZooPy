from setuptools import setup, find_packages

setup(name='pypanda',
    version='0.2',
    description='Python implementation of PUMA, PANDA, LIONESS.',
    url='https://github.com/aless80/pypanda',
    author='Alessandro Marin, Cho-Yi Chen, David van IJzendoorn',
    author_email='AlessandroMarin80@gmail.com',
    license='MIT',
    packages=['pypanda'],
    install_requires=['pandas',
    'numpy',
    'networkx',
    'matplotlib'
    ],
    zip_safe=False)
