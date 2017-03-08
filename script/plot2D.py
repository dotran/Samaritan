from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.pyplot import show, plot




n =int(raw_input())
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_xlim(0,1)
ax.set_ylim(0,1)

while (n!=-1):
    while (n>0):
        xs = float(raw_input())
        ys = float(raw_input())
        n = n-1
        ax.scatter(xs, ys)
    plt.draw()
    plt.show(block=False)   #this creates an empty frozen window.
    n =int(raw_input())
    if(n!=-1):
        plt.clf()
        ax = fig.add_subplot(111)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)


