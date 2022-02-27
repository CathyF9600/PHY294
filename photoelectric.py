from scipy.stats import chisquare,linregress
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
import math as m

def propa_mu1(a,da,f):
    # error propagation for multiplication/division
    df = f*(da/a)
    return df

def propa_mu2(a,da,b,db,f):
    # error propagation for multiplication/division
    df = f*(m.sqrt( (da/a)**2+ (db/b)**2) )
    return df

def propa_avg(L,dL,f):
    # error propagation for taking average
    # L: list of data
    # dL: list of errors
    s = 0
    for i in range(len(L)):
        s += (dL[i]/L[i])**2
    df = f*(m.sqrt(s))
    return df

## Exercise 1

e = 1.60217662 * (10**-19) # electron charge: C
waveL = np.array([935,640,615,590,535,505,455,390]) * (10**-9)
c = 299792458 # speed of light: m/s
f = c/waveL
h = []
b = []
dh = []
db = []

# Metal 1
V1 = np.array([0,0.489,0.526,0.642,0.841,0.957,1.17,1.536]) # stopping volage 1: V
res1 = linregress(f, V1)
slope1 = res1.slope
b.append(res1.intercept)
h.append(slope1 * e) # slope = h/e
print("std",res1.stderr)
dh.append(propa_mu1(res1.slope,res1.stderr,slope1 * e)) # assuming e=const, then dh=h*(dm/m)
db.append(res1.intercept_stderr) # error of intercept

# Metal 2
V2 = np.array([0,0.49,0.525,0.641,0.842,0.956,1.171,1.536])
res2 = linregress(f, V2)
slope2 = res2.slope
b.append(res2.intercept)
h.append(slope2 * e)
dh.append(propa_mu1(res2.slope,res2.stderr,slope2 * e))
db.append(res2.intercept_stderr) # error of intercept

# Metal 3
# V3

# Metal 4
V4 = np.array([0,0.49,0.527,0.642,0.843,0.955,1.17,1.536])
res4 = linregress(f, V4)
slope4 = res4.slope
b.append(res4.intercept)
h.append(slope4 * e)
dh.append(propa_mu1(res4.slope,res4.stderr,slope4 * e))
db.append(res4.intercept_stderr) # error of intercept

# Metal 5
# V5

# Metal 6
# V6

#Find h, f0, E0
h_final = np.average(h)
dh_final = propa_avg(h,dh,h_final)

f0 = -np.divide(np.multiply(b,e),h) # f0 = -be/h, np.multiply do * elementwise, cant use * bc e!=int nor list
df0=[]
for i in range(len(f0)): # cut-off frequency differs from metal to metal
    df0.append(propa_mu2(b[i],db[i],h[i],dh[i],f0[i]))

E0 = h*f0 # same as np.multiply
dE0=[]
for i in range(len(f0)):
    dE0.append(propa_mu2(f0[i],df0[i],h[i],dh[i],E0[i]))

print("Plank's constant h: ", h_final, "±", dh_final)
print("Cut-off frequency: ", f0, "±", df0)
print("Work function: ", E0, "±", dE0)

# graph 

## Exercise 2