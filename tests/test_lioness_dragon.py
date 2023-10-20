import pytest
import os
from netZooPy.lioness.lioness_for_dragon import LionessDragon
import pandas as pd
import numpy as np

def test_lioness_dragon_barcode():
    print("Start LIONESS-DRAGON test by merging on barcode!")
    s = LionessDragon(layer1 = "tests/lioness/dragon/toy_data_layer11.tsv", layer2 = "tests/lioness/dragon/toy_data_layer12.tsv", output_dir = "tests/lioness/dragon/lioness_dragon_results",merge_col = 'barcode',delim="\t")
    s.lioness_loop()
    
def test_lioness_dragon_id():
    print("Start LIONESS-DRAGON test by merging on id!")
    s = LionessDragon(layer1 = "tests/lioness/dragon/toy_data_layer1.tsv", layer2 = "tests/lioness/dragon/toy_data_layer2.tsv", output_dir = "tests/lioness/dragon/lioness_dragon_results",delim="\t")
    s.lioness_loop()
    
def test_lioness_dragon_alldata():
    print("Start LIONESS-DRAGON test with alldata!")
    ext1 = 'expr'
    ext2 = 'meth'
    data1 = pd.read_csv("tests/lioness/dragon/toy_data_layer1.tsv",sep="\t", index_col=0)
    data2 = pd.read_csv("tests/lioness/dragon/toy_data_layer2.tsv",sep="\t", index_col=0)
    data1 = data1.add_suffix(ext1)
    data2 = data2.add_suffix(ext2)
    all_data = pd.merge(data1,data2,left_index=True, right_index=True, how="inner")
    s = LionessDragon(all_data = all_data, ext1=ext1, ext2=ext2, output_dir = "tests/lioness/dragon/lioness_dragon_results_alldata")
    s.lioness_loop()

test_lioness_dragon_barcode()
test_lioness_dragon_id()
test_lioness_dragon_alldata()