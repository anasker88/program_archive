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
    nx = struct.unpack("i", nx_b)[0]
    ny = struct.unpack("i", ny_b)[0]
    u = np.zeros((ny + 2, nx + 1))
    v = np.zeros((ny + 1, nx + 2))
    p = np.zeros((ny + 2, nx + 2))
    t = np.zeros((ny + 2, nx + 2))
    # 以降、時間発展データを読み込む
    while True:
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
            im = plt.imshow(
                t[1 : ny + 1, 1 : nx + 1],
                interpolation="nearest",
                animated=True,
                origin="lower",
                vmax=0.5,
                vmin=-0.5,
            )
            if len(ims) == 0:
                plt.colorbar()
            # ベクトル場を描画
            X = []
            Y = []
            U = []
            V = []
            for i in range(1, nx + 1):
                for j in range(1, ny + 1):
                    X.append(i)
                    Y.append(j)
                    U.append((u[j][i] + u[j][i - 1]) / 2)
                    V.append((v[j][i] + v[j - 1][i]) / 2)
            plt.quiver(X, Y, U, V, angles="xy", scale_units="xy", scale=0.1)
            ims.append([im])
        except EOFError:
            break
        except struct.error:
            break

ani = animation.ArtistAnimation(fig, ims, interval=100)
ani.save("result.gif", writer="pillow")
