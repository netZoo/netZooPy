import pytest
import os

def test_panda():
    python run_panda.py -e ./ToyData/ToyExpressionData.txt -m ./ToyData/ToyMotifData.txt -p ./ToyData/ToyPPIData.txt -o test_panda.txt output_panda.txt