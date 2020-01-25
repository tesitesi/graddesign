# ロケット卒業設計冷却計算

import numpy as np 
import csv 
import matplotlib.pyplot as plt 
import os

from scipy.optimize import curve_fit 


# 全体固定値
"""
# 燃焼ガスについて
Cpg = 6.7224 * 1e+3 #[J/kgK] # スロート部での燃焼ガスの定圧比熱
Cpg *= 0.000947817 / 2.20462 / (9./5.) #[btu/lb degF]

gamma = 1.1294 # [無次元] # 燃焼ガスの比熱比

Pc = 13. #[MPa] # 燃焼圧
Pc *= 145.038 #[psia]

Dt = 0.0422 * 2#[m] # スロート半径
Dt *= 39.3701 #[in]

Rt = 0.2 * Dt #[in] # スロート部曲率半径

cstar = 1770.5 # [m/s] # C*
cstar *= 3.28084 # [fps]

g = 9.80665 #[m/s2] # 重力加速度
g *= 3.28084 #[fps2]

# 冷却剤について
# RP-1
mldot = 10.28 #[kg/s] # 冷却剤流量
mldot *= 2.20462 #[lb/s]

Cpl = 3000. #[J/kgK] # 冷却剤の定圧比熱
Cpl *= 0.000947817 / 2.20462 / (9./5.) #[btu/lb degF]

gammal = 1.2 # 冷却剤の比熱比
"""
# 燃焼ガスについて
Cpg = 7.4257 * 1e+3 #[J/kgK] # スロート部での燃焼ガスの定圧比熱
Cpg *= 0.000947817 / 2.20462 / (9./5.) #[btu/lb degF]

gamma = 1.1259 # [無次元] # 燃焼ガスの比熱比

Pc = 6. #[MPa] # 燃焼圧
Pc *= 145.038 #[psia]

Dt = 0.062 * 2#[m] # スロート直径
Dt *= 39.3701 #[in]

Rt = 0.95 * Dt #[in] # スロート部曲率半径

cstar = 1754.0 # [m/s] # C*
cstar *= 3.28084 # [fps]

g = 9.80665 #[m/s2] # 重力加速度
g *= 3.28084 #[fps2]

# 冷却剤について
# RP-1
mldot = 10.30 #[kg/s] # 冷却剤流量
mldot *= 2.20462 #[lb/s]

# ノズル素材について
#lambdaw = 391 #[W/mK] # 壁内熱伝導係数
lambdaw = 218 / 12. / 3600. #[Btu / F in sec]

e =  1e-3 #[m] 内壁から冷却剤溝までの壁の厚さ
e *= 39.3701 #[in]

a = 4e-3 #[m] #冷却溝高さ
a *= 39.3701 #[in]

b = 4e-3 #[m] #冷却溝幅
b *= 39.3701 #[in]

n = 110.

Al = a*b*n #[m2] # 冷却溝全部の断面積
Al *= 39.3701 * 39.3701 #[in2]


# 計算の初期値
Twg0 = 120. #[K] ：最初の内壁温度
Twg0 *= 9. / 5. #[R]

Tl0 =  100. #[K] :最初の冷却剤温度
Tl0 *= 9. / 5. #[R]

# 出力
out_Tl = []
out_Twg = []


# ==== csvから入力 =====
print('input start')

csv_file = open(os.path.dirname(__file__)+'/cool_input.csv','r',encoding='utf-8-sig')

X = [] # スロートからの距離
R = [] # ノズル半径 
Tg = [] # ガス温度 
M = [] # マッハ数
P = [] # 燃焼ガス圧力
for row in csv.reader(csv_file):
    X.append(row[0])
    R.append(row[1])
    Tg.append(row[2])
    M.append(row[3])
    P.append(row[4])
X = np.array(X, dtype=np.float64)[::-1] #[m]
R = np.array(R,dtype=np.float64)[::-1] #[m]
Tg = np.array(Tg,dtype=np.float64)[::-1] #[K]
M = np.array(M,dtype=np.float64)[::-1] #[無次元]
P = np.array(P,dtype=np.float64)[::-1] #[MPa]

R *= 39.3701 #[in]
Tg *= 9. / 5. #[R]
P *= 145.038 #[psia]


print('input complete')
# ==== csv読み込みここまで =======




# ==== 熱伝達係数の計算式を定義 ====
print('h definition start')

# RP-1 の粘性係数を温度（ランキン）に対して返す。
def u(T):
    #if T <= 500:
        #u = 1.743 / 0.41337887 /  3600. / 12.
    #else:
    #u =  4.91664334e-07  * np.exp(2.78880905e+03 / T)
    u = 0.5 / 0.41337887 /12. / 3600.# [lb / in s]
    return u


# RP-1 の動粘性係数を in2/s で返す。
def nyu(T):
    #if T <= 500:
    #    nyu = 0.7678 * 0.0393701 * 0.0393701 * 5
    #else:
    #nyu = 2.77251253e-01 / (T - 4.46325281e+02)
    nyu = 0.0033573067 #[in2/s]
    return nyu 

# RP-1 の熱伝導率を返す。
def lambdal(T):
    #lambdal = 2.11068937e-03 / (T + 7.69870886e+02) -9.35918399e-08
    lambdal = 0.05781759824 / 12. / 3600. #[btu / degF in s]
    return lambdal

# RP-1の定圧比熱を実験値を２次式に近似してフィッティング。
def Cpl(T):
    #Cpl = 3.78841696e-07 * T * T + (-3.47542440e-05) * T + 3.97452633e-01
    Cpl = 1.194229483135 #[btu / lb degF]
    return Cpl
"""
# RP-1の熱伝導率を近似式で返す。
def kl(T):
    kl =  -9.21389611e-10  * T +  2.07981274e-06

    return kl
"""
# ガス側熱伝達係数計算 #MKS単位系であることに注意
def hg(Tg, Twg, M, R):
    ug = 46.6e-10 * (26.6**0.5) * (Tg**0.6)
    Tcns = Tg * (1. + 0.5 * (gamma - 1.) * M * M)
    Prg = 4. * gamma / (9. * gamma - 5.)
    Pcns = Pc
    sigma = (0.5 * Twg / Tcns * (1. + 0.5 * (gamma - 1.) * M * M) + 0.5)**(-0.68) * (1. + 0.5 * (gamma - 1.) * M * M) **(-0.12)
    At = np.pi * Dt * Dt / 4.
    A = np.pi * R * R
    print(R/(Dt/2.))
    hg = 0.026 / Dt**0.2 *(ug**0.2 * Cpg / Prg**0.6) * (Pcns* g / cstar)**0.8 * (Dt / Rt)**0.1 * (At / A)**0.9 * sigma
    return hg

# ヌッセルト数と冷却剤の熱伝導率から熱伝達率を計算する
def hc(Tl, Twl):
    Prl = u(Tl) * Cpl(Tl) / lambdal(Tl)
    d = a * np.sqrt(2) #[in] #水力直径
    
    A_coolpath = Al / n #[in] #冷却溝一個の断面積
    rho_loc = 0.0252891 #RP-1の密度, [lb/in3]
    Re = (mldot/n) * d / rho_loc / A_coolpath / nyu(Tl)
    Nu = 0.0214 * (Re**0.8) * (Prl**0.4) * (u(Tl) / u(Twl))**0.14 # Sieder-Tateの式
    hc = Nu * lambdal(Tl) / d
    hc *= n #冷却管の数だけ足す
    """
    G = mldot / Al #[lb/in2sec]
    hc = 0.029 * Cpl(Tl) * u(Tl)**0.2 / Prl**(2./3.) * (G**0.8 / d**0.2) * (Tl/Twl)**0.55
    hc *= n
    """
    return hc 


print('h definition complete')
# ===========================  熱伝達係数の計算式ここまで ===================================




# =========================== ループ計算ここから ===================================
print('Loop Calc. Start')
# 初期値
Tl = Tl0
Twg = Twg0
# 刻み幅ごとに進行していくループ計算
for x in range(len(R)):
    if np.abs(X[x]) <= 0.06:
        maxitr = 500
    else:
        maxitr = 500
    if x % 50 == 0:
        print('Looping at x = ' + str(x))
    if x != 0:
        Tl = Tl_next
        Twg = Twg_next
    # ２つの方法で計算した Tl_next が一致するまで計算する
    epsilon = 1000.
    count = 0
    while epsilon > 1. and count < (2 * maxitr) :
        Twg = Twg + 1.
        # 収束するまで1Kずつずらす
        # 燃焼ガスからノズル内壁への熱流束
        q = hg(Tg[x], Twg, M[x], R[x]) * (Tg[x] - Twg)
        # x+1との冷却剤温度差
        A_danmen = 0.0393701 * 2. * R[x] * np.pi #[in2]
        deltaTl = q * A_danmen / (mldot *Cpl(Tl))
        Tl_next = Tl + deltaTl
        # x+1の壁面温度をxの内壁温度と同じに設定
        Twg_next = Twg
        # 壁内熱伝導の式から、x+1での外壁温度を出す
        Twl_next = Twg_next - q / lambdaw * e
        # 冷却剤と外壁感の熱伝達の式から、 x+1 での冷却剤温度を出す
        Tl_next_2 = Twl_next - q / hc(Tl_next, Twl_next)
        count = count + 1
        epsilon = np.abs(Tl_next - Tl_next_2)
        if count == (2 * maxitr):
            print('Warning: Loop Broke Out at X_from_throat = ' + str(X[x]) + " [m]!")
    # データ保存
    out_Tl.append(Tl)
    out_Twg.append(Twg)


print('Loop Calc. Complete!')
# ============================== ループ計算ここまで =======================================-



    
# ======================= 出力ここから =====================
print('csv output start')

out_Twg = np.array(out_Twg) #[R]
out_Tl = np.array(out_Tl) #[R]
out_Twg *= 5. / 9. #[K]
out_Tl *= 5. / 9.#[K]
np.savetxt(os.path.dirname(__file__)+'/cool_out_Twg.csv', out_Twg.T)
np.savetxt(os.path.dirname(__file__)+'/cool_out_Tl.csv', out_Tl.T)

print('csv output complete')
# ==================== 出力ここまで ==================--

# =================== グラフ描画ここから ===============-
print('graphing start')

plt.title("Temperature in Nozzle")
plt.xlabel('Distance from Throat [m]')
plt.ylabel('Temperature [K]')
plt.scatter(X,out_Tl, label='Tl',s=1)
#plt.scatter(X,out_Twg, label='Twg',s=1)
plt.legend()
plt.show()

# =================== グラフ描画ここまで ==================

print('End of Program')



        



        