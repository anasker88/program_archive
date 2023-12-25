import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

np.set_printoptions(threshold=np.inf)

# Pokemon.csvから種族値データを読み込み、dataに格納する
data = []
with open("Pokemon.csv", "r") as f:
    for line in f:
        data.append(line.strip().split(","))
# dataを整形
# 1行目はヘッダーなので削除
data.pop(0)
# 種族値のみを抽出
data = [d[5:11] for d in data]
# 数値に変換
original_data = [[int(d) for d in l] for l in data]
# 平均を0にする
original_data = np.array(original_data)
data = original_data - np.mean(original_data, axis=0)
# dataをPCAで二次元に
# 　散布行列を求める
C = data.T.dot(data)
print("散布行列", C)
# 　固有値と固有ベクトルを求める
eig_val, eig_vec = sp.linalg.eig(C)
# 固有値の大きい順にならべかえる
idx = np.argsort(eig_val)[::-1]
eig_val = eig_val[idx]
eig_vec = eig_vec[:, idx]
# 上から2つの固有ベクトルを取り出す
eig_vec = eig_vec[:, :2]
print("pcaの固有ベクトル\n", eig_vec)
# データを変換する
data_dim2 = np.array(data)
data_dim2 = np.dot(data, eig_vec)
# k-means法でクラスタリングする
# 3つのクラスタに分ける
# 　初期値をランダムに設定
# 　データの最小値と最大値を求める
min_x = np.min(data_dim2[:, 0])
max_x = np.max(data_dim2[:, 0])
min_y = np.min(data_dim2[:, 1])
max_y = np.max(data_dim2[:, 1])
# 　初期値をランダムに設定
c1 = (
    np.random.rand() * np.array([max_x - min_x, 0])
    + np.random.rand() * np.array([0, max_y - min_y])
    + np.array([min_x, min_y])
)
c2 = (
    np.random.rand() * np.array([max_x - min_x, 0])
    + np.random.rand() * np.array([0, max_y - min_y])
    + np.array([min_x, min_y])
)
c3 = (
    np.random.rand() * np.array([max_x - min_x, 0])
    + np.random.rand() * np.array([0, max_y - min_y])
    + np.array([min_x, min_y])
)

# 　クラスタリングを行う
for i in range(10):
    cluster1 = []
    cluster2 = []
    cluster3 = []
    for d in data_dim2:
        # 　距離を計算
        dist1 = np.linalg.norm(d - c1)
        dist2 = np.linalg.norm(d - c2)
        dist3 = np.linalg.norm(d - c3)
        # 　最も近いクラスタに追加
        if dist1 < dist2 and dist1 < dist3:
            cluster1.append(d)
        elif dist2 < dist1 and dist2 < dist3:
            cluster2.append(d)
        else:
            cluster3.append(d)
    # 　クラスタの重心を計算
    c1 = np.mean(cluster1, axis=0)
    c2 = np.mean(cluster2, axis=0)
    c3 = np.mean(cluster3, axis=0)
# 　クラスタを描画
import matplotlib.pyplot as plt

cluster1 = np.array(cluster1)
cluster2 = np.array(cluster2)
cluster3 = np.array(cluster3)
# クラスタを描画
plt.scatter(cluster1[:, 0], cluster1[:, 1], c="red")
plt.scatter(cluster2[:, 0], cluster2[:, 1], c="blue")
plt.scatter(cluster3[:, 0], cluster3[:, 1], c="green")
# 　重心を描画
plt.scatter(c1[0], c1[1], c="red", marker="x")
plt.scatter(c2[0], c2[1], c="blue", marker="x")
plt.scatter(c3[0], c3[1], c="green", marker="x")
# pca_and_kmeans.pngという名前で保存
plt.savefig("pca_and_kmeans.png")


# kernel PCA
# 　カーネル関数を定義(ガウシアンカーネル)
def kernel(x, y):
    return np.exp(-np.linalg.norm(x - y) ** 2 / 10**6)


# 　カーネル行列を計算
K = np.zeros((len(original_data.T), len(original_data.T)))
# print(original_data)
for i in range(6):
    for j in range(6):
        K[i, j] = kernel(original_data.T[i], original_data.T[j])
# 中心化
H = np.identity(len(original_data.T)) - 1 / len(original_data) * np.ones(
    (len(original_data.T), len(original_data.T))
)
K = np.dot(np.dot(H, K), H)
print("カーネル行列\n", K)
#  PCAと同様に固有値と固有ベクトルを求める
eig_val, eig_vec = sp.linalg.eig(K)
# 固有値の大きい順にならべかえる
idx = np.argsort(eig_val)[::-1]
eig_val = eig_val[idx]
eig_vec = eig_vec[:, idx]
# 上から3つの固有ベクトルを取り出す
eig_vec = eig_vec[:, :3]
print("kernel pcaの固有ベクトル\n", eig_vec)
# データを変換する
data_dim3 = np.array(original_data)
data_dim3 = np.dot(original_data, eig_vec)
# print("kernel pcaのデータ", data_dim3)
# k-means法でクラスタリングする
min_x = np.min(data_dim3[:, 0])
max_x = np.max(data_dim3[:, 0])
min_y = np.min(data_dim3[:, 1])
max_y = np.max(data_dim3[:, 1])
min_z = np.min(data_dim3[:, 2])
max_z = np.max(data_dim3[:, 2])
# 適当にデータを3つのクラスタに分ける
cluster1 = data_dim3[: len(data_dim3) // 3]
cluster2 = data_dim3[len(data_dim3) // 3 : len(data_dim3) // 3 * 2]
cluster3 = data_dim3[len(data_dim3) // 3 * 2 :]


# 　クラスタリングを行う
def dist(x, cluster):
    n = len(cluster)
    return (
        kernel(x, x)
        - 2 / (n + 0.001) * np.sum([kernel(x, y) for y in cluster])
        + 1 / (n + 0.001) ** 2 * np.sum([kernel(y, y) for y in cluster])
    )


print("cluster1", len(cluster1))
print("cluster2", len(cluster2))
print("cluster3", len(cluster3))
for i in range(10):
    tmp1 = []
    tmp2 = []
    tmp3 = []
    for d in data_dim3:
        # 　距離を計算
        dist1 = dist(d, cluster1)
        dist2 = dist(d, cluster2)
        dist3 = dist(d, cluster3)
        # 　最も近いクラスタに追加
        if dist1 < dist2 and dist1 < dist3:
            tmp1.append(d)
        elif dist2 < dist1 and dist2 < dist3:
            tmp2.append(d)
        else:
            tmp3.append(d)
    # 　クラスタを更新
    cluster1 = np.array(tmp1)
    cluster2 = np.array(tmp2)
    cluster3 = np.array(tmp3)
    print("step", i, "done")
    print("cluster1", len(cluster1))
    print("cluster2", len(cluster2))
    print("cluster3", len(cluster3))

# 　クラスタを描画
plt.figure()
if len(cluster1) != 0:
    c1 = np.mean(cluster1, axis=0)
    plt.scatter(cluster1[:, 0], cluster1[:, 1], c="red")
    plt.scatter(c1[0], c1[1], c="red", marker="x")
if len(cluster2) != 0:
    c2 = np.mean(cluster2, axis=0)
    plt.scatter(cluster2[:, 0], cluster2[:, 1], c="blue")
    plt.scatter(c2[0], c2[1], c="blue", marker="x")
if len(cluster3) != 0:
    c3 = np.mean(cluster3, axis=0)
    plt.scatter(cluster3[:, 0], cluster3[:, 1], c="green")
    plt.scatter(c3[0], c3[1], c="green", marker="x")
# 　重心を描画
# kernel_pca_and_kmeans.pngという名前で保存
plt.savefig("kernel_pca_and_kmeans.png")
