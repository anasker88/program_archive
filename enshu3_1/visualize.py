import struct

import matplotlib.animation as animation
import matplotlib.pyplot as plt

fig = plt.figure()
ims = []


# データファイルを開く
with open("result", "rb") as f:
    # 最初はnx,nyを読み込む
    nx_b = f.read(4)
    ny_b = f.read(4)
    nx = struct.unpack("i", nx_b)[0]
    ny = struct.unpack("i", ny_b)[0]
    u = [[0 for i in range(nx + 1)] for j in range(ny + 2)]
    v = [[0 for i in range(nx + 2)] for j in range(ny + 1)]
    p = [[0 for i in range(nx + 2)] for j in range(ny + 2)]
    t = [[0 for i in range(nx + 2)] for j in range(ny + 2)]
    # 以降、時間発展データを読み込む
    while True:
        try:
            # 最初はcycleを読み込む
            cycle = f.read(4)
            # uを読み込む
            for i in range(nx + 1):
                for j in range(ny + 2):
                    data = f.read(8)
                    u[j][i] = struct.unpack("d", data)
            # vを読み込む
            for i in range(nx + 2):
                for j in range(ny + 1):
                    data = f.read(8)
                    v[j][i] = struct.unpack("d", data)
            # pを読み込む
            for i in range(nx + 2):
                for j in range(ny + 2):
                    data = f.read(8)
                    p[j][i] = struct.unpack("d", data)
            # tを読み込む
            for i in range(nx + 2):
                for j in range(ny + 2):
                    data = f.read(8)
                    t[j][i] = struct.unpack("d", data)
            # 　温度は色で表現
            #  端のデータは仮想セルのものなので表示しない
            im = plt.imshow(
                t[1 : nx + 1][1 : ny + 1],
                interpolation="nearest",
                animated=True,
            )
            # ベクトル場を描画
            # X = []
            # Y = []
            # U = []
            # V = []
            # for i in range(1, nx + 1):
            #     for j in range(1, ny + 1):
            #         X.append(i)
            #         Y.append(j)
            #         U.append(u[j][i][0])
            #         V.append(v[j][i][0])
            # plt.quiver(X, Y, U, V, angles="xy", scale_units="xy", scale=0.1)
            ims.append([im])
        except EOFError:
            break
        except struct.error:
            break

ani = animation.ArtistAnimation(fig, ims, interval=100)
ani.save("result.gif", writer="pillow")
