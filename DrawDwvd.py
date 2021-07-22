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
    resolution = int(firstLine[1])

    lines = lines[1:]

    matrix = []
    for line in lines :
        arr = line.split(" ")
        vals = []
        for val in arr :
            if val == "" :
                continue
            vals.append(float(val))

        matrix.append(vals)

    matdf = pd.DataFrame(matrix)
    plt.figure()
    sns.heatmap(matdf)
    plt.savefig(args[2])
