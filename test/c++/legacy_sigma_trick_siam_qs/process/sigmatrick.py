#!/usr/bin/env python2
# Python script for computing the spectral function using the self-energy
# trick
# Rok Zitko, zitko@theorie.physik.uni-goettingen.de, June 2008
# QSZ-version, rok.zitko@ijs.si, March 2009

import os
import math

pi = math.pi

# Spectral functions
Afn = "Ad.dat"
Ffn = "Bd.dat"

# Kramer's Kronig transforms
reAfn = 're' + Afn
reFfn = 're' + Ffn

def kk(fnin, fnout) : os.system('kk ' + fnin + ' ' + fnout)

kk(Afn, reAfn)
kk(Ffn, reFfn)

# SIAM parameters
delta = 0.0
U = 1.0
Gamma = 0.06
epsilond = delta-U/2
B = 0.0
epsilond = epsilond - B/2

def first(x)  : return float( x.split() [0] )
def second(x) : return float( x.split() [1] )

# Energy mesh
f = open(Afn, 'r')
l = f.readlines()
xvals = map(first,l)
f.close()

# Hybridisation function defintion.
# Here it is assumed to be indpedent of energy.
def iif(cond, trueval, falseval):
    if cond:
        return trueval
    else:
        return falseval

yvals = [iif(abs(x)<=1.0, Gamma, 0.0) for x in xvals]

# Write the hybdrisation function.
f = open('imdos.dat', 'w')
for i in range(len(xvals)) :
  str = '%18.15f %18.15f\n' % (xvals[i] * 1.0, yvals[i])
  f.write(str)
f.close()

# Perform the Kramer-Kronig transformation
fnim = 'imdos.dat'
fnre = 'redos.dat'
kk(fnim, fnre)

# Read the results
def readsecond(fn) :
    f = open(fn, 'r')
    l = f.readlines()
    vals = map(second, l)
    f.close()
    return vals

imGammavals = readsecond(fnim)
reGammavals = readsecond(fnre)
os.unlink(fnim)
os.unlink(fnre)

imAvals = readsecond(Afn)
reAvals = readsecond(reAfn)
os.unlink(reAfn)

imFvals = readsecond(Ffn)
reFvals = readsecond(reFfn)
os.unlink(reFfn)

# Calculate the spectral function
f = open('spec.dat', 'w')

Gammavals = []
for i in range(len(xvals)):
  omega = xvals[i]
  # Raw Green's function
  Gold = 1j*imAvals[i] + reAvals[i]
  # F function
  F = 1j*imFvals[i] + reFvals[i]
  # Self-energy
  Sigma = U * (F / Gold + 0.5)
  # (complex) hybridisation
  Gamma = 1j*imGammavals[i] + reGammavals[i]
  # Green's function
  G = 1.0/(omega - epsilond + Gamma - Sigma)
  # Spectral density
  A = -1.0/pi * G.imag
  str = '%18.10f %18.15f\n' % (xvals[i], A)
  f.write(str)
  
f.close()
