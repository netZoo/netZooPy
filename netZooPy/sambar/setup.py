from setuptools import setup, find_packages

setup(name='pysambar',
    version='0.1',
    description='Python implementation of SAMBAR',
    url='https://github.com/genisott/pysambar',
    author='Genís Calderer',
    author_email='genis.calderer@gmail.com',
    license='MIT',
    packages=['pysambar'],
    install_requires=['pandas',
    'numpy',
    'networkx',
    'matplotlib',
    'scipy'
    ],
    package_data={'pysambar': ['ToyData/*']},
    #data_files=[('ToyData', ["ToyData/esizef.csv","ToyData/genes.txt","ToyData/mut.ucec.csv","ToyData/h.all.v6.1.symbols.gmt"])],
    zip_safe=False)
