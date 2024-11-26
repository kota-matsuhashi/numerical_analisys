import matplotlib.pyplot as plt

# データの準備（各iterごとにリストに格納）
data = {
    0: [1.0, 0.4605877, 0.1262056, 0.0691241, 0.0345620, 0.0184062, 0.0092031, 0.0048511, 0.0024256, 0.0012757],
    1: [1.0, 0.4346645, 0.1701164, 0.0891916, 0.0403540, 0.0241375, 0.0117835, 0.0065551, 0.0034194, 0.0017995],
    2: [1.0, 0.4249923, 0.1934927, 0.0854295, 0.0349410, 0.0220533, 0.0095037, 0.0058501, 0.0026726, 0.0015968],
    3: [1.0, 0.4665783, 0.1989545, 0.0893916, 0.0388877, 0.0203000, 0.0095915, 0.0054348, 0.0027868, 0.0015898],
    4: [1.0, 0.4588463, 0.1950968, 0.0908113, 0.0379258, 0.0197206, 0.0090472, 0.0049990, 0.0024928, 0.0014434],
    5: [1.0, 0.4628320, 0.1905254, 0.0863033, 0.0363492, 0.0196543, 0.0088828, 0.0049154, 0.0029830, 0.0013966],
    6: [1.0, 0.4691027, 0.1881280, 0.0871967, 0.0379901, 0.0190368, 0.0118658, 0.0055250, 0.0039949, 0.0018669],
    7: [1.0, 0.4728940, 0.1911197, 0.0876227, 0.0466144, 0.0210679, 0.0146588, 0.0068149, 0.0049217, 0.0022965],
    8: [1.0, 0.4704120, 0.2025849, 0.0907654, 0.0533913, 0.0243340, 0.0169892, 0.0078884, 0.0056923, 0.0026529],
    9: [1.0, 0.4729076, 0.2208333, 0.0987129, 0.0582313, 0.0269020, 0.0187870, 0.0087144, 0.0062842, 0.0029260]
}

# 各 na に対応するインデックス
na_values = list(range(1, 11))

# グラフのプロット
plt.figure(figsize=(10, 6))

for iter_num, rmax_values in data.items():
    plt.plot(na_values, rmax_values, label=f'Iter {iter_num+1}')

# グラフの装飾
plt.xlabel('na')
plt.ylabel('rmax')
plt.title('rmax values for each iteration')
plt.yscale('log')  # rmaxの範囲が広いので対数スケールを使用
plt.legend()
plt.grid(True)

# グラフの表示
plt.show()