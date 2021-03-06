# 卒業設計
# ベルノズルのコンター（半径分布）を決める
import numpy as np 
import matplotlib.pyplot as plt 

X = []
contour = []
dx = 0.001

r_t = 0.062 # スロート部半径
r_e = 0.6515 # 出口部半径
l_conb = 0.193 # 燃焼室長さ
theta_i = 30.0 # 膨張接続部角度
theta_e = 5.0 # 出口角度

theta_i = theta_i * np.pi / 180.
theta_e = theta_e * np.pi / 180.

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


#b2 = (l_d * np.tan(theta_i) - x1 * np.tan(theta_e)) / (np.tan(theta_i) - np.tan(theta_e))
#a2 = (np.tan(theta_i)-np.tan(theta_e)) / (2*l_d - 2*x1)

a2 = (np.tan(theta_i)-np.tan(theta_e)) / (2*l_d - 2*x1)
b2 = np.tan(theta_e) + 2 * a2 * l_d


for x in [round(i * 0.001, 3) for i in range(int((x1)*1000), int((l_d+dx)*1000), int(dx*1000))]:
    #contour.append(-a2 * (x - b2)**2 + a2 * (l_d - b2)**2 + r_e)
    if -a2 * x * x + b2 * x + r_t > r_e :
        break
    X.append(x)
    contour.append(-a2 * x * x + b2 * x + r_t)
    


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

#np.savetxt('nozzlecontour_output.csv', contour.T)
np.savetxt('nozzlecontour_output.csv', X.T)
