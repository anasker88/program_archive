import struct

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ims = []


# データファイルを開く
with open("result", "rb") as f:
    # 最初はnx,nyを読み込む
    nx_b = f.read(4)
    ny_b = f.read(4)
    data_num_b = f.read(4)
    nx = struct.unpack("i", nx_b)[0]
    ny = struct.unpack("i", ny_b)[0]
    data_num = struct.unpack("i", data_num_b)[0]
    u = np.zeros((ny + 2, nx + 1))
    v = np.zeros((ny + 1, nx + 2))
    p = np.zeros((ny + 2, nx + 2))
    t = np.zeros((ny + 2, nx + 2))

    # 以降、時間発展データを読み込む
    def updatefig(frame):
        # プログレスバー
        print(f"\rReading data... {frame+1}/{data_num}", end="")
        # 全てクリア
        plt.clf()
        try:
            # 最初はcycleを読み込む
            cycle = f.read(4)
            # uを読み込む
            for i in range(nx + 1):
                for j in range(ny + 2):
                    data = f.read(8)
                    (u[j][i],) = struct.unpack("d", data)
            # vを読み込む
            for i in range(nx + 2):
                for j in range(ny + 1):
                    data = f.read(8)
                    (v[j][i],) = struct.unpack("d", data)
            # pを読み込む
            for i in range(nx + 2):
                for j in range(ny + 2):
                    data = f.read(8)
                    (p[j][i],) = struct.unpack("d", data)
            # tを読み込む
            for i in range(nx + 2):
                for j in range(ny + 2):
                    data = f.read(8)
                    (t[j][i],) = struct.unpack("d", data)
            # 　温度は色で表現
            #  端のデータは仮想セルのものなので表示しない
            # y軸は上下逆になっているので注意
            plt.imshow(
                t[1 : ny + 1, 1 : nx + 1],
                interpolation="nearest",
                animated=True,
                origin="lower",
                vmax=0.5,
                vmin=-0.5,
            )
            plt.colorbar()
            # ベクトル場を描画
            X = []
            Y = []
            U = []
            V = []
            W = []
            for i in range(1, nx + 1):
                for j in range(1, ny + 1):
                    X.append(i - 1)
                    Y.append(j - 1)
                    a = (u[j][i] + u[j][i - 1]) / 2
                    b = (v[j][i] + v[j - 1][i]) / 2
                    w = 0.0
                    if abs(a) > abs(b):
                        w = a * np.sqrt(1 + (b / a) ** 2)
                    elif b != 0:
                        w = b * np.sqrt(1 + (a / b) ** 2)
                    else:
                        w = 0
                    if w > 0:
                        U.append(a / w)
                        V.append(b / w)
                        W.append(w)
                    else:
                        U.append(0)
                        V.append(0)
                        W.append(0)
            # 色の濃さで速度を表現
            plt.quiver(
                X, Y, U, V, W, cmap="Greys", scale=nx * 1.5, pivot="mid", angles="xy"
            )
        except EOFError:
            return
        except struct.error:
            return

    ani = animation.FuncAnimation(fig, updatefig, interval=100, frames=data_num)
    ani.save("result.gif", writer="pillow")
