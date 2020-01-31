# 卒業設計
# ベルノズル半径分布をもとに、各部の流れマッハ数を算出する
import numpy as np 
import matplotlib.pyplot as plt 
import csv

csv_file = open('/Users/tesiyosi/dev/graddesign/nozzle/nozzlemach_input.csv','r',encoding='utf-8-sig')

X = [] 
A = []
for row in csv.reader(csv_file):
    X.append(row[0])
    A.append(row[1])
X = np.array(X,dtype=np.float64)
A = np.array(A,dtype=np.float64)

# CEAの結果
M_inj = 0.
M_combend = 0.407
M_throat = 1.
M_exit = 4.345

gamma_inj = 1.1344
gamma_combend = 1.1332
gamma_throat = 1.1294
gamma_exit = 1.1823

gamma = 0.

# dx によって変わる
x_inj = 0
x_combend =  277
x_throat = 341
x_exit = 1079

M = np.ndarray(x_exit+1)

def mach(A_prev,A, M, gamma):
    dA = A - A_prev
    left = 2 * dA * (1 + 0.5 * (gamma - 1) * M * M) / (( M * M - 1) * (A + A_prev))
    dM = M * left
    return M + dM


# 燃焼室部
for x in range(0, x_combend, 1):
    M[x] = M_combend / x_combend * x

# 収縮部
for x in range(x_combend, x_throat, 1):
    #gamma = (gamma_combend * (x-x_throat) + gamma_throat * (x_combend-x)) / (x_combend - x_throat)
    #M[x] = mach(A[x-1],A[x],M[x-1],gamma)
    M[x] = (M_combend * (x_throat - x) + M_throat * (x - x_combend)) / (x_throat - x_combend)

M[x_throat] = 1. + 1e-2

# 膨張部
for x in range(x_throat+1, x_exit, 1):
    gamma = (gamma_throat * (x-x_exit) + gamma_exit * (x_throat-x)) / (x_throat - x_exit)
    M[x] = mach(A[x-1],A[x],M[x-1],gamma)



plt.title("Mach Number in Nozzle")
plt.xlabel('Distance from Throat [m]')
plt.ylabel('Mach Number [1]')
print(X.shape)
print(M.shape)
plt.scatter(X,M,s=2,c='red')
plt.grid(True)
plt.show()

np.savetxt('/Users/tesiyosi/dev/graddesign/nozzle/nozzlemach_output.csv', M.T)
