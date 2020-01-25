import numpy as np 
import matplotlib.pyplot as plt 
import csv

csv_file = open('nozzletemp_input.csv','r',encoding='utf-8-sig')

X = [] 
M = []
for row in csv.reader(csv_file):
    X.append(row[0])
    M.append(row[1])
X = np.array(X,dtype=np.float64)
M = np.array(M,dtype=np.float64)

# CEAの結果
M_inj = 0.
M_combend = 0.407
M_throat = 1.
M_exit = 3.41795

gamma_inj = 1.1309
gamma_combend = 1.1297
gamma_throat = 1.1259
gamma_exit = 1.1747

T_inj = 3678.34
T_combend = 3636.42
T_throat = 3497.34
T_exit = 1939.54


P_inj = 6.000
P_combend = 5.0538
P_throat = 3.2139
P_exit = 0.004674

gamma = 0.


x_inj = 1
x_combend = 192
x_throat = 286
x_exit = 2046

T_0 = T_throat * (1 + 0.5 * (gamma_throat - 1 ))
P_0 = P_throat * (1 + 0.5 * (gamma_throat - 1 )) **(gamma_throat/(gamma_throat-1))


T = []
P = []


# 燃焼室部
for x in range(0, x_combend, 1):
    gamma = (gamma_inj * (x - x_inj) + gamma_combend * (x_combend - x)) / (x_combend - x_inj)
    print(gamma)
    T.append( T_0 / (1 + 0.5 * (gamma - 1 )* M[x]*M[x]))
    P.append(P_0 / (1 + 0.5 * (gamma - 1 )* M[x]*M[x])**(gamma/(gamma-1)))

# 収縮部
for x in range(x_combend, x_throat, 1):
    gamma = (gamma_combend * (x-x_throat) + gamma_throat * (x_combend-x)) / (x_combend - x_throat)
    T.append( T_0 / (1 + 0.5 * (gamma - 1 )* M[x]*M[x]))
    P.append(P_0 / (1 + 0.5 * (gamma - 1 )* M[x]*M[x])**(gamma/(gamma-1)))

# 膨張部
for x in range(x_throat, x_exit+1, 1):
    gamma = (gamma_throat * (x-x_exit) + gamma_exit * (x_throat-x)) / (x_throat - x_exit)
    T.append( T_0 / (1 + 0.5 * (gamma - 1 )* M[x]*M[x]))
    P.append(P_0 / (1 + 0.5 * (gamma - 1 )* M[x]*M[x])**(gamma/(gamma-1)))
"""
plt.title("Gas Temperature in Nozzle")
plt.xlabel('Distance from Throat [m]')
plt.ylabel('Temperature [K]')
plt.scatter(X,T,s=2,c='green')
plt.grid(True)
plt.show()
"""
plt.title("Gas Pressure in Nozzle")
plt.xlabel('Distance from Throat [m]')
plt.ylabel('Pressure [MPa]')
plt.scatter(X,P,s=2,c='blue')
plt.grid(True)
plt.show()

P = np.array(P)
T = np.array(T)

np.savetxt('nozzletemp_output.csv', T.T)