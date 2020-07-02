import csv 
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt 
def eq(y, t, N, beta, gamma, delta, alpha, rho):
    S, E, I, R, D = y 
    dSdt = -beta(t) * S * I / N 
    dEdt = beta(t) * S * I / N - delta * E 
    dIdt = delta * E - (1 - alpha) * gamma * I - alpha * rho * I 
    dRdt = (1 - alpha) * gamma * I 
    dDdt = alpha * rho * I 
    return dSdt, dEdt, dIdt, dRdt, dDdt
N = 330000000

D = 10.5 #Infection Period
gamma = 1.0 / D 
delta = 1.0 / 4.0  # Incubation Time
R_0s = 5.0 
k = .0345
x0 = 50
R_0f = 1.0 
def DynamicR_0(t):
    return (R_0s-R_0f) / (1 + np.exp(-k*(-t+x0))) + R_0f
def beta(t):
    return DynamicR_0(t) * gamma
alpha = 0.005  # death rate
rho = 1/5  # Infection -> Death
S0, E0, I0, R0, D0 = N-40, 20, 20, 0, 0  # initial conditions
t = np.linspace(0, 160, 160) 
y0 = S0, E0, I0, R0, D0

# Integrate the SIR equations over the time grid, t.
ret = odeint(eq, y0, t, args=(N, beta, gamma, delta, alpha, rho))
S, E, I, R, D = ret.T
I1 = [] # Actual Infected
D1 = [] # Actual Dead
R1 = [] # Actual Recovered
with open('daily.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
    for row in spamreader:
        a = row[0]
        a = a.split(",")
        I1.append(a[2])
        D1.append(a[13])
        R1.append(a[11])
D1.pop(0)
R1.pop(0)
I1.pop(0)
I1 = np.array(I1)
I1 = I1.astype(np.float)
I1 = I1[::-1]
D1 = np.array(D1)
D1[D1 == ''] = 0 
D1 = D1.astype(np.float)
D1 = D1[::-1]
R1 = np.array(R1)
R1[R1 == ''] = 0 
R1 = R1.astype(np.float)
R1 = R1[::-1]
F1 = I1 - R1 - D1
I = I[30::]
F1 = F1[30::]
t2 = np.linspace(0, 130, 130)
def plotseird(t, S, E, I, R, D): 
  #plt.plot(t, S, 'b', alpha=0.7, linewidth=2, label='Susceptible')
  #plt.plot(t, E, 'y', alpha=0.7, linewidth=2, label='Exposed')
  plt.plot(t2, I, 'r', alpha=0.7, linewidth=2, label='Predicted I(t)')
  #plt.plot(t, R, 'g', alpha=0.7, linewidth=2, label='Recovered')
  #plt.plot(t, D, 'k', alpha=0.7, linewidth=2, label='Dead')
  plt.plot(t2, F1, 'o', alpha=0.7, linewidth=2, label='Actual I(t)')
  plt.xlabel('Days Since March 1')
  plt.legend()
  plt.title("COVID-19 SIRD/SEIR Model Infected Compartment")
  plt.ylim(0,2000000)
  plt.show()
plotseird(t, S, E, I, R, D)
