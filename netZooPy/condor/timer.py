"""A tic-toc analog in Python
Adapted from http://stackoverflow.com/questions/5849800/tic-toc-functions-analog-in-python
"""
import time

class Timer(object):
    def __init__(self, name=None,silent=False):
        self.silent = silent
        if self.silent: return
        if name:
            print(name)
        

    def __enter__(self):
        self.tic = time.time()

    def __exit__(self, type, value, traceback):
        if self.silent: return
        print('  Elapsed time: %.2f sec.' % (time.time() - self.tic))
