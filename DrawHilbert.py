import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

X = []
Y = []
with open("hilbertReal.txt", "r") as f:
    lines = f.read().splitlines()
    for line in lines:
        xy = line.split(" ")
        X.append(float(xy[0]))
        Y.append(float(xy[1]))

#X, Y = zip(*sorted(zip(X,Y)))

plt.plot(X,Y)
plt.savefig('hilbertReal.png')