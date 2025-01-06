import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ファイルからデータを読み込む
data = np.loadtxt('velocity.DAT', skiprows=1)

# データを分割する
u = data[:, 0].reshape((21, 21, 21))
v = data[:, 1].reshape((21, 21, 21))
w = data[:, 2].reshape((21, 21, 21))

# 格子点の座標を生成する
x = np.linspace(0, 2, 21)
y = np.linspace(0, 2, 21)
z = np.linspace(0, 2, 21)
X, Y, Z = np.meshgrid(x, y, z)

# 速度ベクトルのプロット
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.quiver(X, Y, Z, u, v, w, length=0.1)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Velocity Field')
plt.show()