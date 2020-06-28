import numpy as np
import matplotlib.pyplot as plt 
import csv 
import pandas as pd
import requests
import io
import urllib.request
import random
## Start of Initial Setup with United States' Data
url = 'https://raw.githubusercontent.com/nytimes/covid-19-data/master/us.csv'
read_data=requests.get(url).content
address=pd.read_csv(io.StringIO(read_data.decode('utf-8')))
address = np.array(address)
url1 = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv'
read_data1=requests.get(url1).content
address1=pd.read_csv(io.StringIO(read_data1.decode('utf-8')))
address1 = np.array(address1)
Date = address[:,0][-1]
N = 330000000
I0 = address[:,1][-1]
R0 = address1[225][-1]
D0 = address[:,2][-1]
S0 = N - I0 - R0
t = np.linspace(0, 200) 
y0 = S0,I0, R0, D0
beta = 0.10      #infectivity rate
gamma = 1/21    #recovery rate
delta = .0049   #mortality rate
## End of Initial Setup
def equations(y, t, N, beta, gamma, delta):
    S, I, R, D = y 
    dSdt = -beta * S * I / N 
    dIdt = beta * S * I / N - gamma * I - delta * I 
    dRdt = gamma * I 
    dDdt = delta * I 
    return dSdt, dIdt, dRdt, dDdt
#Integrate
fin = odeint(equations, y0, t, args=(N, beta, gamma, delta))
S, I, R, D = fin.T

def plotseird(t, S, I, R, D): 
  plt.plot(t, S, 'chartreuse', linewidth = 3, label='Susceptible')
  plt.plot(t, I, 'red', linewidth = 3, label='Infected')
  plt.plot(t, R, 'blue', linewidth = 3, label='Recovered')
  plt.plot(t, D, 'orchid', linewidth = 3, label='Dead')
  plt.xlabel('Days after {Date}'.format(Date = Date))
  plt.ylabel('Count')
  plt.legend()
  plt.title('SIRD Model of COVID-19 for U.S.')
  plt.show()
plotseird(t, S, I, R, D)
