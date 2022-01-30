from scipy.stats import chisquare,linregress
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker

x=np.array([2.45E-01,
2.24E-01,
2.07E-01,
1.94E-01,
1.83E-01,
1.73E-01,
1.65E-01,
1.58E-01]) # lamda: angstrom

# radius ring1
r1_in = np.array([13.15, 12.45, 11.55, 10.25, 9.5, 9.3, 8.9, 8.75])
r1_out = np.array([15.95, 14.8, 12.35, 11.25, 10.65, 11.1, 10.55, 10.5])
y1 = (r1_in+r1_out)/2
print(y1)

# radius ring2
r2_in = np.array([22.75, 21.8, 18.95, 17.25, 16.45, 16.3, 15.95, 15.15])
r2_out = np.array([27.1, 25.15, 21.75, 21.25, 20.35, 19.3, 18.25, 17.55])
y2 = (r2_in+r2_out)/2

# linear fit ring1
res1 = linregress(x, y1)
print(linregress(x, y1))
R=66 #mm
d1 = 2*R*(1/res1.slope)
print("m1:", res1.slope,"±",res1.stderr)
print(f"R-squared: {res1.rvalue**2:.6f}")

# linear fit ring2
res2 = linregress(x, y2)
print(linregress(x, y2))
d2 = 2*R*(1/res2.slope)
print("m2:", res2.slope,"±",res2.stderr)
print(f"R-squared: {res2.rvalue**2:.6f}")

# linear lines
y1e=res1.intercept + res1.slope*x
y2e=res2.intercept + res2.slope*x

# error bars
xerr_bars = len(y1)*[0.005] # angstrom
yerr_bars = len(y1)*[0.01767766953] # mm

## start plotting!
plt.xlabel('wavelength(angstrom)')
plt.ylabel('radius(mm)')
plt.title('Curve Fitting')
# axis limit
plt.xlim(0.14, 0.26)
plt.ylim(7.75, 30)

# plot data
plt.scatter(x,y1,s=5,c='k',label='$r=61\lambda$')
plt.errorbar(x,y1,yerr=yerr_bars,xerr=xerr_bars,fmt='None',ecolor='b')
# plot linear fit
plt.plot(x,y1e,'--r',label='Linear Fit for d_100')

# plot data
plt.scatter(x,y2,s=5,c='k',label='$r=100\lambda$')
plt.errorbar(x,y2,yerr=yerr_bars,xerr=xerr_bars,fmt='None',ecolor='r')
# plot linear fit
plt.plot(x,y2e,'--b',label='Linear Fit for d_110')

xticks = np.arange(0.16,0.26,0.02)
yticks = np.arange(10,32.5,2.5)

plt.legend()
plt.show()

# print("chi1:", chisquare(y1, y1e, ddof=6, axis=0).statistic)
# print("chi2:", chisquare(y2, y2e, ddof=6, axis=0).statistic)

# find reduced chi
def chisqr(obs, exp, error):
    chisqr = 0
    for i in range(len(obs)):
        chisqr = chisqr + ((obs[i]-exp[i])**2)/(error[i]**2)
    return chisqr/6

print("y1-axis chi sq:", chisqr(y1, y1e, yerr_bars))
print("y2-axis chi sq:", chisqr(y2, y2e, yerr_bars))
