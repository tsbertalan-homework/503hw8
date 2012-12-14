import numpy as np
exp = np.exp
ln = np.log
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(10, 5))
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)

R = 8.3145
NA = 6.022e-23
kb = R / NA
e = 120 * kb
s = 0.34e-3
T = 300
resolution = 10000
rmin = 3.32e-4
rmax = 8.0e-4
dr = (rmax - rmin) / resolution

g = lambda r: exp(-4 * e / kb / T * ((s/r) ** 12 - (s/r) ** 6))
zeroThis = lambda r: -4 * e / kb / T * (-12*s**12/r**13 + 6*s**6/r**7)
gp = lambda r: zeroThis(r) * g(r)


rl = list(np.arange(rmin, rmax, dr))
gl = [g(r) for r in rl]
gpl = [gp(r) for r in rl]
CDFold = 0

ax1.plot(rl, gl, 'k')
ax2.plot(rl, gpl, 'k')
fig.suptitle('Hw. 8, No. 1', fontsize=16)
ax1.set_title(r'pair correlation function $g(r)$')
ax1.set_xlabel(r'Separation $[m]$')
ax2.set_xlabel(r'Separation $[m]$')

#ax.set_ylabel(r'Potential $[J]$')

ymin = min(gl)
ymax = max(gl)

absyl = [abs(y) for y in gpl]
zeroIndex = absyl.index(min(absyl))
gpl0 = rl[zeroIndex]
ax2.set_title(r'$\partial g(r) / \partial r$')
ax2.axhline(y=0, color='k')

ax2.scatter([gpl0], [0], color='k')
ax1.scatter([gpl0], max(gl), color='k')
ax2.annotate(r'$\partial g(r) / \partial r = 0$ at $r=%.3f [nm]$'%(gpl0*1000), xy=(gpl0+.1*gpl0, .1*max(gpl)))

ax1.set_xlim((min(rl), max(rl)))
ax2.set_xlim((min(rl), max(rl)))
#ax1.set_ylim((ymin - (ymax - ymin)/10, ymax))
#ax2.set_xlim((gpl0, max(rl)))
#ax2.set_ylim((min(gpl), max(gpl)))

print rl[0]
print ymin
print rl[-1]
#ax1.legend()
fig.tight_layout()

from os.path import basename
fig.savefig(basename(__file__) + '.pdf')
fig.savefig(basename(__file__) + '.png')
plt.show()
