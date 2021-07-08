import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

X = []
Y = []
with open("analyzed.txt", "r") as f:
    lines = f.read().splitlines()
    for line in lines:
        yx = line.split(" ")
        Y.append(float(yx[0]))
        X.append(float(yx[1]))

X, Y = zip(*sorted(zip(X,Y)))

plt.plot(X,Y)
plt.savefig('fft.png')