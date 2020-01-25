from scipy.optimize import curve_fit 
import numpy as np 
import csv 
import matplotlib.pyplot as plt 
def f(x,a,b,c):
    return a / (x + b) + c

def f2(x,a,b):
    return a / (x + b)


def f3(x, a,b,c):
    return a * x * x + b * x + c

def f4(x,a,b):
    return a * np.exp(b / x)

 
# RP-1 の粘性係数を温度（ランキン）に対して返す。
def u(T):
    t_sample = np.array([373.15, 363.15, 353.15, 343.15, 333.15, 323.18, 313.12, 303.18, 293.38]) #[K]
    t_sample = t_sample * 9. / 5.
    u_sample = np.array([0.5727, 0.6362, 0.7103, 0.8007, 0.9100, 1.050, 1.225, 1.447, 1.743]) #[mPa s]
    u_sample = u_sample / 0.41337887 /  3600. / 12. # [lb/in s]
    opt, _ = curve_fit(f4, t_sample, u_sample)
    print(opt)
    u =  opt[0] * np.exp(opt[1] / T)
    return u

# RP-1 の動粘性係数を in2/s で返す。
def nyu(T):
    t_sample = np.array([373.15, 363.15, 353.15, 343.15, 333.15, 323.18, 313.12, 303.18, 293.38]) #[K]
    t_sample = t_sample * 9. / 5.
    nyu_sample = np.array([0.7678, 0.8445, 0.9335, 1.042, 1.173, 1.341, 1.549, 1.814, 2.166]) #[mm2/s]
    nyu_sample = nyu_sample * 0.0393701 * 0.0393701 #[in2/s]
    opt, _ = curve_fit(f2, t_sample, nyu_sample)
    nyu =  opt[0] / (T + opt[1])
    print(opt)
    return nyu 

# RP-1 の熱伝導係数を返す。
def lambdal(T):
    # 10 MPa
    # https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=901440
    t_sample = np.array([293.15, 300.65, 317.15, 333.19, 339.17, 353.20, 353.17, 359.15, 374.67,
    393.68, 421.85, 453.05, 475.83, 500.89, 533.05, 533.18, 591.26, 613.23, 693.08, 732.68]) #[K]
    t_sample = t_sample * 9. / 5.
    sigmal_sample = np.array([0.114, 0.115, 0.110, 0.107, 0.108, 0.106, 0.106, 0.103, 0.102, 0.100, 0.096, 
    0.093, 0.090, 0.089, 0.085, 0.081, 0.081, 0.076, 0.072, 0.068]) #[W/mK]
    sigmal_sample = sigmal_sample * 0.0478778 / 3600 # [(Btu/lb F )* lb/in sec] = [Btu / F in sec]
    opt, _ = curve_fit(f, t_sample, sigmal_sample)
    print(opt)
    sigmal =  opt[0] / (T + opt[1]) + opt[2]
    return sigmal

# RP-1の定圧比熱の実験値
def Cpl(T):
    # 10MPa
    t_sample = np.array([293.76, 334.15, 373.42, 434.65, 475.45, 535.32, 576.63, 633.84, 671.42]) #[K]
    t_sample = t_sample * 9. / 5. #[R]
    Cpl_sample = np.array([2015., 2153., 2305., 2531., 2699., 2965., 3229., 3565., 3810.]) #[J/kgK]
    Cpl_sample *= 0.000947817 / 2.20462 / (9./5.) #[btu/lb degF]
    opt, _ = curve_fit(f3, t_sample, Cpl_sample)
    print(opt)
    Cpl = opt[0] * T**2 + opt[1] * T + opt[2]
    return Cpl


def k(T):
    # https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=902380
    t_sample = np.array([300, 350, 400, 450]) #[K]
    t_sample = t_sample * 9. / 5. #[R]
    k_sample = np.array([0.118, 0.113, 0.105, 0.10]) # [W/mK]
    k_sample *= 0.04815 / 3600 #[btu/ in s degF]
    opt, _ = curve_fit(f3, t_sample, k_sample)
    print(opt)
    k = opt[0] * T**2 + opt[1] * T + opt[2]
    return k



t_sample = np.array([373.15, 363.15, 353.15, 343.15, 333.15, 323.18, 313.12, 303.18, 293.38]) #[K]
t_sample = t_sample * 9. / 5.
nyu_sample = np.array([0.7678, 0.8445, 0.9335, 1.042, 1.173, 1.341, 1.549, 1.814, 2.166]) #[mm2/s]
nyu_sample = nyu_sample * 0.0393701 * 0.0393701 #[in2/s]
T = np.arange(500, 5000)
plt.plot(T, nyu(T))
plt.scatter(t_sample,nyu_sample)
plt.show()
