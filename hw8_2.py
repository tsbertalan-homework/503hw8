# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 20:28:04 2012

@author: tsbertalan
"""
import numpy as np
exp = np.exp
ln = np.log
pi = np.pi
import matplotlib.pyplot as plt

ax = range(4)
fig = range(4)

fig[0] = plt.figure(figsize=(10, 8))
ax[0] = fig[0].add_subplot(1, 1, 1)

fig[1] = plt.figure(figsize=(10, 8))
ax[1] = fig[1].add_subplot(1, 1, 1)

fig[2] = plt.figure(figsize=(10, 8))
ax[2] = fig[2].add_subplot(1, 1, 1)

resolution = 100
dt = (5. - 1) / resolution
tl = list(np.arange(1, 5, dt))
def Eandg(ksi, l):
    
    z = 20
    
    def eta_f(gl):
        eta_ = 0
        for (t, g) in zip(tl, gl):
            eta_ += 4 * pi * t ** 2 * g * dt / z
        return eta_
        
    def Qin_f(r, gl):
        print "Qin_f"
        Qin_ = z
        count = 0
        for (t, g) in zip(tl, gl):
            if t <= r:
                Qin_ -= 4 * pi * t ** 2 * g * dt
                cout += 1
        print count
        return Qin_
    
    def Ueff_f(r, gl):
        sum_ = -z * ksi
        for (t, g) in zip (tl, gl):
            if t >= r:
                sum_ += 4 * pi * ksi * t * g * dt
        return sum_
        
    def gl_f(r, gl):
        return exp(-Ueff_f(r, gl)) / eta_f(gl)
        
    gl = [20./50. for x in range(resolution)]
    gold = gl
    threshold = 1e-5
    normtest = 1
    while normtest > threshold:
        normtest = abs(np.linalg.norm(np.array(gl) - np.array(gold)))*10000
        gstar = [gl_f(r, gl) for r in tl]
        gl = [l * gs * (1-l) * go for (gs, go) in zip(gstar, gold)]
        gold = gl
    def E(gl):
        sum_ = 0
        for (t, g) in zip(tl, gl):
            sum_ += t * g * dt
        return -z * 4 * pi * sum_

    return (E(gl), gl)
    
for (ksi, l, i) in zip([1/5., 1./50, 1./2], [0.1, .5, 0.02], range(4)):
    (energy, pairCorrelation) = Eandg(ksi, l)
    pairCorrelation.reverse()
    ax[i].plot(tl, pairCorrelation)
    ax[i].set_title(r'$\xi = %.2f$'%ksi)
    ax[i].set_xlabel(r'r (units of $\sigma$)')
    ax[i].set_ylabel(r'$g(r)$')
    print "for ksi =", ksi, ", energy =", energy
    filename = 'hw8_2_%i'%i
    fig[i].savefig(filename + '.png')
    fig[i].savefig(filename + '.pdf')
plt.show()