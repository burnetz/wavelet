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
    N = int(firstLine[0])
    maxLevel = int(firstLine[1])

    lines = lines[1:]
    i = 0
    inputData = []
    while i < maxLevel :
        inputData.append([[],[]])
        #print(inputData)
        #print(i)
        i += 1
    
    fig = plt.figure(figsize=(6.0,10.0))   
    i = 0 
    maxVal = 0
    while i < N - N/pow(2, maxLevel):
        #print(i)
        lxy = lines[i].split(" ")
        l = int(lxy[0]) - 1
        inputData[l][0].append(float(lxy[1]))
        inputData[l][1].append(float(lxy[2]))
        maxVal = max(maxVal, abs(float(lxy[2])))
        #print(inputData[l])
        i += 1
    
    #print(inputData)
    i  = 0
    while i < maxLevel :
        #print(i)
        ax = fig.add_subplot(maxLevel, 1, i + 1)
        ax.stem(inputData[i][0], inputData[i][1], use_line_collection=True)
        ax.set_xlim(0,N)
        ax.set_ylim(-maxVal*1.1,maxVal*1.1)
        i += 1

    plt.savefig(args[2])
