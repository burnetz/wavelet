# 参考 : https://qiita.com/futurebone/items/fd42892c237a817cdbb4

import numpy as np
import wave
import struct
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ファイルを読み出し
#print("Please Enter the filename")
#wavf = input("->")
wavf = "./asset/voice.wav"
wr = wave.open(wavf, 'r')

# waveファイルが持つ性質を取得
ch = wr.getnchannels()
width = wr.getsampwidth()
fr = wr.getframerate()
fn = wr.getnframes()

print("Channel: ", ch)
print("Sample width: ", width)
print("Frame Rate: ", fr)
print("Frame num: ", fn)
print("Params: ", wr.getparams())
print("Total time: ", 1.0 * fn / fr)

# waveの実データを取得し、数値化
data = wr.readframes(wr.getnframes())
wr.close()

# フレーム数は2のべき乗にする
fn = 512
X = np.arange(0,fn,1)
Y = np.frombuffer(data, dtype=np.int16)
Y1 = Y[0:fn]/(2.0**(width*8))

# 読み込み結果をグラフで確認したいとき
# plt.plot(X,Y1)
# plt.savefig('figure.png')

# テキストファイルに出力
with open('wave.txt', 'w') as f:
    f.write("%d\n" % fn)
    for d in Y1:
        f.write("%.12f\n" % d)