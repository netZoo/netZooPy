## Description
Forked from QuackenbushLab/pypanda.
I work on run_puma and pypanda/puma.py. Those methods can still have pandas instead of puma, same for this README.  

## TODO
Rewrite this README. Anaconda nice to mention. Update the example from pandas to puma. update citations and descriptions. check all links. figures!  

Check David's python commands in this README  

## PyPuma (Python Puma)
Python implementation of PUMA (PANDA Using MicroRNA Associations)  

PANDA: 
_Glass K, Huttenhower C, Quackenbush J, Yuan GC. Passing Messages Between Biological Networks to Refine Predicted Interactions, PLoS One, 2013 May 31;8(5):e64832_

### Table of Contents
* [Panda implementation](#panda-algorithm)  
* [Installation](#installation)  
* [Usage](#usage)  
* [python](#run-from-python)
* [Terminal](#run-from-the-terminal)
* [Results](#results)

### Panda algorithm
To find agreement between the three input networks first the responsibility (R) is calculated.  

<img src="https://github.com/aless80/pypanda/raw/develop/img/responsibility.png" height="30">  

Thereafter availability (A) is calculated.  

<img src="https://github.com/aless80/pypanda/raw/develop/img/availability.png" height="30">  

Availability and responsibility are combined with the following formula.  

<img src="https://github.com/aless80/pypanda/raw/develop/img/combine.png" height="30">  

Protein cooperativity and gene co-regulatory networks are updated.  

<img src="https://github.com/aless80/pypanda/raw/develop/img/cooperativity.png" height="30">  
<img src="https://github.com/aless80/pypanda/raw/develop/img/co-regulatory.png" height="30">  

P and C are updated to satisfy convergence.  

<img src="https://github.com/aless80/pypanda/raw/develop/img/p.png" height="30">  
<img src="https://github.com/aless80/pypanda/raw/develop/img/c.png" height="30">  

Hamming distance is calculated every iteration.  

<img src="https://github.com/aless80/pypanda/raw/develop/img/hamming.png" height="40">  


### Installation
PyPuma requires Python 2.7. We recommand the following commands to install PyPuma (on Ubuntu and Debian derived systems, also works on OSX):
#### Using a virtual environment
Using [python virtual environment](http://docs.python-guide.org/en/latest/dev/virtualenvs/) is the cleanest installation method. 

Cloning git and setting up the [python virtual environment](http://docs.python-guide.org/en/latest/dev/virtualenvs/):
```no-highlight
pip install --user pipenv   #Make sure you have Pipenv
git clone https://github.com/aless80/pypanda.git
cd pypuma
virtualenv pypumaenv #Create a folder for the virtual environment inside the cloned git folder 
source pypumaenv/bin/activate
```
Installing pypuma:
```no-highlight
(pypumaenv)$ pip install -r requirements.txt
(pypumaenv)$ python setup.py install #--user
```

Complete uninstall of pypuma:
```no-highlight
(pypuma)$ deactivate	#Exit
rm -rf pypumaenv
write about uninstalling the setup using --record
```

#### Using pip 
Using pip on the user's install directory
```no-highlight
git clone https://github.com/aless80/pypanda.git
cd pypuma
python setup.py install --user
#to run from the command line you will need to make pypuma executable and add the bin directory to your PATH:
cd bin
chmod +x pypuma
echo "$(pwd):PATH" >> ~/.bashrc
source ~/.bashrc
```
To run PyPuma from Windows (tested on Windows 10) install Git (https://git-scm.com/downloads) and Anaconda Python2.7 (https://www.continuum.io/downloads) and from the Anaconda Prompt run:
```no-highlight
git clone https://github.com/aless80/pypanda.git
cd pypuma
python setup.py install
```
### Usage
#### Run from the terminal
PyPandas can be run directly from the terminal with the following options:
```
-h help
-e, --expression (required) expression values
-m, --motif (optional) pair file of motif edges, when not provided analysis continues with Pearson correlation matrix
-p, --ppi (optional) pair file of PPI edges
-f, --remove (optional) remove missing values (default is Fales)
-o, --output (required) output file
-q, --outputlioness (optional) output lioness single sample network
-s, --s1 (required) s1
-x, --s2 (required) s2
-t, --t1 (required) t1
-z, --t2 (required) t2
```
To run PyPuma on the example data:
```
pypanda -e ToyData/ToyExpressionData.txt -m ToyData/ToyMotifData.txt -p ToyData/ToyPPIData.txt -f True -o test_puma.txt -q test_lioness.txt -s 0 -x 1 -t 1 -z 2
```
To reconstruct a single sample Lioness Pearson correlation network:
```
$ pypanda -e ToyData/ToyExpressionData.txt -o test_puma_pearson.txt -q test_lioness_pearson.txt
```
#### Run from python
Fire up your python shell or ipython notebook. 
Import PyPuma library:
```python
from pypuma import Puma
from pypuma import Lioness
import pandas as pd
```
Run Puma algorithm, leave out motif and PPI data to use Pearson correlation network:
```python
p = Puma('ToyData/ToyExpressionData.txt', 'ToyData/ToyMotifData.txt', 'ToyData/ToyPPIData.txt', remove_missing=False)
```
Save the results:
```python
p.save_puma_results(file = 'Toy_Puma.pairs')
```
Return a network plot:
```python
plot = AnalyzePuma(p)
plot.top_network_plot(top=100, file='top_100_genes.png')
```
Calculate indegrees for further analysis:
```python
indegree = p.return_puma_indegree()
```
Calculate outdegrees for further analysis:
```python
outdegree = p.return_puma_outdegree()
```
Run the Lioness algorithm for single sample networks:
```python
l = Lioness(p)
```
Save Lioness results:
```python
l.save_lioness_results(file = 'Toy_Lioness.txt')
```
Return a network plot for one of the Lioness single sample networks:
```python
plot = AnalyzeLioness(l)
plot.top_network_plot(column= 0, top=100, file='top_100_genes.png')
```
### Results
```
Example Panda output:
TF  Gene  Motif Force
---------------------
CEBPA	AACSL	0.0	-0.951416589143
CREB1	AACSL	0.0	-0.904241609324
DDIT3	AACSL	0.0	-0.956471642313
E2F1	AACSL	1.0	3.6853160511
EGR1	AACSL	0.0	-0.695698519643

Example lioness output:
Sample1 Sample2 Sample3 Sample4
-------------------------------
-0.667452814003	-1.70433776179	-0.158129613892	-0.655795512803
-0.843366539284	-0.733709815256	-0.84849895139	-0.915217389738
3.23445386464	2.68888472802	3.35809757371	3.05297381396
2.39500370135	1.84608635425	2.80179804094	2.67540878165
-0.117475863987	0.494923925853	0.0518448588965	-0.0584810456421

TF, Gene and Motif order is identical to the panda output file.
```
