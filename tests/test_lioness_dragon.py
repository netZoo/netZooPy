import pytest
import os
from netZooPy.lioness.lioness_for_dragon import LionessDragon
import pandas as pd
import numpy as np

def test_lioness_dragon():
    print("Start LIONESS-DRAGON test!")

s = LionessDragon(layer1 = "lioness/dragon/toy_data_layer1.tsv", layer2 = "lioness/dragon/toy_data_layer2.tsv", output_file = "lioness/dragon/lioness_dragon_results.csv")
s.lioness_loop()
