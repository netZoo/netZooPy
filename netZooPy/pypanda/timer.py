"""A tic-toc analog in Python

Adapted from http://stackoverflow.com/questions/5849800/tic-toc-functions-analog-in-python
"""
import time

class Timer(object):
    def __init__(self, name=None):
        if name:
            print(name)

    def __enter__(self):
        self.tic = time.time()

    def __exit__(self, type, value, traceback):
        print('  Elapsed time: %.2f sec.' % (time.time() - self.tic))
