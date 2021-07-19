import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

X1 = []
Y1 = []
with open("hilbertReal.txt", "r") as f:
    lines = f.read().splitlines()
    for line in lines:
        xy = line.split(" ")
        X1.append(float(xy[0]))
        Y1.append(float(xy[1]))

#X, Y = zip(*sorted(zip(X,Y)))

plt.title("Real")
plt.plot(X1,Y1)
plt.savefig('hilbertReal.png')

#plt.close()
X2 = []
Y2 = []
with open("hilbertImag.txt", "r") as f:
    lines = f.read().splitlines()
    for line in lines:
        xy = line.split(" ")
        X2.append(float(xy[0]))
        Y2.append(float(xy[1]))

plt.title("Real, Imag")
plt.plot(X2,Y2)
plt.savefig('hilbertImag.png')

X3 = []
Y3 = []
with open("hilbertPower.txt", "r") as f:
    lines = f.read().splitlines()
    for line in lines:
        xy = line.split(" ")
        X3.append(float(xy[0]))
        Y3.append(float(xy[1]))

plt.title("Real, Imag, power")
plt.plot(X3,Y3)
plt.savefig('hilbertPower.png')

X4 = []
Y4 = []
with open("hilbertFreq.txt", "r") as f:
    lines = f.read().splitlines()
    for line in lines:
        xy = line.split(" ")
        X4.append(float(xy[0]))
        Y4.append(float(xy[1]))

plt.title("Real, Imag, power, freq")
plt.plot(X4,Y4)
plt.savefig('hilbertFreq.png')