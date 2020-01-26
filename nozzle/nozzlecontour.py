# ロケット卒業設計
# ベルノズルのコンター（半径分布）を決める

import numpy as np 
import matplotlib.pyplot as plt 

X = [] # 結果を入れる行列
contour = [] #結果を入れる行列

dx = 0.001 # 格子間隔 [m]

r_t = 0.062 # スロート部半径 [m]
epsilon = 110. #開口比

r_e = r_t * np.sqrt(epsilon) # 出口部半径 [m]
l_conb = 0.193 # 燃焼室長さ [m]
theta_i = 30.0 # 膨張接続部角度
#theta_e = 8.0 # 出口角度 #廃止

theta_i = theta_i * np.pi / 180.
#theta_e = theta_e * np.pi / 180.

r_conb = np.sqrt(3.0) * r_t # 燃焼室＝ノズル入口半径

rk_c = 1.5 * r_t # 収縮部曲率半径
rk_d = 0.4 * r_t # 膨張部曲率半径

# まずは収縮部
a_c = 1.0 / (2.0 * rk_c) # 収縮部二次関数の係数

l_c = np.round(np.sqrt((r_conb - r_t) / a_c),3) # 収縮部長さ

# 燃焼室部半径
for x in [round(i * 0.001, 3) for i in range(int((-l_c-l_conb)*1000), int(-l_c*1000), int(dx*1000))]:
    X.append(x)
    contour.append(r_conb)

# 収縮部半径
for x in [round(i * 0.001, 3) for i in range(int(-l_c*1000), int(0*1000), int(dx*1000))]:
    X.append(x)
    contour.append(a_c * x * x + r_t)

# つづいて膨張部
a_d1 = 1.0 / (2.0 * rk_d) # スロート直後の二次関数（下に凸）
l_d = round(0.80 * (r_e - r_t) / np.tan(15. * np.pi / 180.),3) # 膨張部全体長さ

x1 = np.round(np.tan(theta_i) / (2 * a_d1),2) # 膨張接続部座標
# スロート部
X.append(0.00)
contour.append(r_t)
# スロート直後
for x in [round(i * 0.001, 3) for i in range(int(0+dx)*1000, int(x1*1000), int(dx*1000))]:
    X.append(x)
    contour.append(a_d1 * x * x + r_t)


# theta_i, theta_e を指定して、出口径を無視する場合
#a2 = - (np.tan(theta_i)-np.tan(theta_e)) / (2*l_d - 2*x1)
#b2 = np.tan(theta_e) + 2 * a2 * l_d

# 新しい：出口径、theta_iを指定する。
x2 = x1 + (r_e - r_t) / np.tan(15 * np.pi / 180.) #ノズル出口座標
y1 = a_d1 * x1 * x1 + r_t

a2 = ((x1 - x2) + (r_e - y1)) / (x1 - x2)**2
b2 = 0.9 - 2 * a2 * x1
c2 = y1 + a2 * x1 * x1 - x1


#for x in [round(i * 0.001, 3) for i in range(int((x1)*1000), int((l_d+dx)*1000), int(dx*1000))]:
for x in [round(i * 0.001, 3) for i in range(int((x1)*1000), int((x2)*1000), int(dx*1000))]:
    X.append(x)
    contour.append(a2 * x * x + b2 * x + c2)
    if (a2 * x * x + b2 * x + c2 > r_e):
        break
    
theta_e = np.arctan((contour[-1] - contour[-2]) / dx) * 180. / np.pi
print('最後の角度：' + str(theta_e))

X = np.array(X)
contour = np.array(contour)

plt.title("Nozzle Contour")
plt.axes().set_aspect('equal', 'datalim')
plt.xlabel('Length [m]')
plt.ylabel('Nozzle Radius [m]')
plt.scatter(X,contour,s=2, c='blue')
plt.scatter(X,-contour,s=2,c='blue')
plt.grid(True)
plt.show()

np.savetxt('/Users/tesiyosi/dev/graddesign/nozzle/nozzlecontour_output_contour.csv', contour.T)
np.savetxt('/Users/tesiyosi/dev/graddesign/nozzle/nozzlecontour_output_X.csv', X.T)
