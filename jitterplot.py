import os
import re
import sys
import matplotlib.pyplot as plt
import numpy as np

def rand_jitter(arr):
    stdev = .03*(max(arr)-min(arr))
    return arr + np.random.randn(len(arr)) * stdev

def jitter(x, y, c, s=0.07, marker='o'):
   # plt.figure()
    return plt.scatter(rand_jitter(x), rand_jitter(y), s=s, c=c, marker=marker)


