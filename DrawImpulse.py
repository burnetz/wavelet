import sys
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

args = sys.argv

with open(args[1], "r") as f:
    lines = f.read().splitlines()
    firstLine = lines[0].split(" ")
    maxLevel = int(firstLine[0])
    N = int(firstLine[1])

    lines = lines[1:]
    i = 0
    inputData = []

    while i < maxLevel:
        inputData.append([[],[]])
        #print(inputData)
        #print(i)
        i += 1
    
    fig = plt.figure(figsize=(10.0,4.0))   
    i = 0 
    maxVal = 0
    while i < N * maxLevel :
        #print(i)
        lxy = lines[i].split(" ")
        l = i // N
        #if abs(float(lxy[1])) > 0.00001 :
        inputData[l][0].append(float(lxy[0]))
        inputData[l][1].append(float(lxy[1]))
        #print(inputData[l])
        i += 1
    
    #print(inputData)
    i  = 0
    while i < maxLevel:
        #print(i)
        ax = fig.add_subplot(maxLevel, 1, i + 1)
        ax.stem(inputData[i][0], inputData[i][1], use_line_collection=True, markerfmt=' ')
        ax.set_xlim(0,N)
        ax.set_ylim(-max(inputData[i][1])*1.1, max(inputData[i][1])*1.1)
        i += 1

    plt.savefig(args[2])
