"""
variance of wrong model parameters 
and parameters correlation
"""
import random
from math import*
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import*

def test_func(x, v0, dv, a1, a2, a3, b1, b2):
    return (a1+a2*(x-v0)+a3*sin(x-v0))/((x-v0)**2+dv**2)+b1*(x-v0)+b2*(x-v0)**2

def fit_func(x, v0, dv,a,b,c,d,f):
    return (a+b*(x-v0)+c*(x-v0)**2)/((x-v0)**2+dv**2)+d*(x-v0)+f*(x-v0)**2

# model parameters
v0 = 6.114 # center freq
dv = 0.5  # spectral width
N = 100  # array length
kNum = 500  # statistics for fitting
level = 0.5  # noise level
a1 = 1.2;
a2 = 0.5;
a3 = 1.5;
b1 = 0.03;
b2 = 0.005;
p0 = (v0, dv, a1, a2, a3, b1, b2)
   
x = np.linspace(v0 - 10 * dv,v0 + 10 * dv, N)
test = np.linspace(0, 0, N)
fit = np.linspace(0, 0, N)
for i in range(N):
    test[i] = test_func(x[i], v0, dv, a1, a2, a3, b1, b2)
    fit[i] = fit_func(x[i], v0, dv, a1, a2, a3, b1, b2)
test_f = np.linspace(0,0,N)
fit_f = np.linspace(0,0,N)
mid_w = np.linspace(0,0,7)
#diff = np.zeros((kNum, 7))
# correlation matrix for wrong model parameters
corr = np.zeros((7, 7)) 
matr_w = np.zeros((kNum,7))
matr = np.zeros((kNum,7))
# parameters variance 
s_w = np.zeros((7, kNum)) # for wrong model
s = np.zeros((7, kNum)) # for fitting model
for i in range(kNum): 
    for j in range(N):
        test_f[j] = test[j] + random.random()*level
        fit_f[j] = fit[j] + random.random()*level
    p = curve_fit(fit_func, x, test_f)[0]   #[v1,dv1,a,b,c,d,f]=curve_fit(...)[0]
    p1 = curve_fit(fit_func, x, fit_f)[0]
    p[1]=abs(p[1])
    p1[1]=abs(p1[1])
    for k in range(7):
        mid_w[k] = mid_w[k] + p[k]
    for k in range(7):
        #diff[i][k]=abs(p0[k]-m)
        matr_w[i][k] = p[k]
        matr[i][k] = p1[k]
        sum_w = 0
        sum_ = 0
        for l in range(i):
            sum_w = sum_w + (matr_w[l][k] - p0[k])**2
            sum_ = sum_ + (matr[l][k] - p0[k])**2
        sum_w = sum_w / (i + 1)
        s_w[k][i] = sum_w
        sum_ = sum_ / (i + 1)
        s[k][i] = sum_   
mid_w = mid_w / kNum
fout = open("wr_corr.txt","w")
for i in range(7):
    fout.write('\n')
    for j in range(7):
        s0 = 0
        if i == j:
            corr[j][i] = 1
            fout.write('1.000'+' ')
        if i != j:
            for k in range(kNum):
                s0 += (matr_w[k][i] - mid_w[i]) * (matr_w[k][j] - mid_w[j])
            corr[j][i] = s0 / kNum / sqrt(s_w[i][kNum - 1]) / sqrt(s_w[j][kNum - 1])
            fout.write(str("%.3f"%(corr[j][i]))+' ')
fout.close()

par_names=['v0','dv','a1','a2','a3','b1','b2']
n = np.linspace(1, kNum, kNum)
for i in range(7):
    fig = plt.figure(figsize = (7.65, 5.25))
    plt.plot(n, s_w[i], label = u'неправильная модель', color = 'r')
    plt.plot(n, s[i], label = u'правильная модель', color = 'b')
    plt.xlabel('N');
    plt.ylabel('дисперсия ' + par_names[i]);
    plt.legend()
