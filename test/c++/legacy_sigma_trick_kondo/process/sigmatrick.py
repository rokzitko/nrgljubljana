#!/usr/bin/env python2
# Python script for computing the spectral function using the self-energy
# trick
# Rok Zitko, zitko@theorie.physik.uni-goettingen.de, June 2008

import os
import math

pi = math.pi

# Spectral functions
Afn = "Af.dat"
Ffn = "OAf.dat"

# Kramer's Kronig transforms
reAfn = 're' + Afn
reFfn = 're' + Ffn

def kk(fnin, fnout) : os.system('kk ' + fnin + ' ' + fnout)

kk(Afn, reAfn)
kk(Ffn, reFfn)

# Kondo parameters
J = 0.2

def first(x)  : return float( x.split() [0] )
def second(x) : return float( x.split() [1] )

# Energy mesh
f = open(Afn, 'r')
l = f.readlines()
l = filter(lambda x : x.split() [0] != "#", l)
xvals = map(first, l)
f.close()

# Hybridisation function defintion.
# Here it is assumed to be indpedent of energy.
def iif(cond, trueval, falseval):
    if cond:
        return trueval
    else:
        return falseval

yvals = [iif(abs(x)<=1.0, 0.5, 0.0) for x in xvals]

# Write rho.
f = open('imrho.dat', 'w')
for i in range(len(xvals)) :
    str = '%22.18f %18.15f\n' % (xvals[i] * 1.0, yvals[i])
    f.write(str)

f.close()
        
# Perform the Kramer-Kronig transformation
fnim = 'imrho.dat'
fnre = 'rerho.dat'
kk(fnim, fnre)

# Read the results
def readsecond(fn) :
    f = open(fn, 'r')
    l = f.readlines()
    l = filter(lambda x : x.split() [0] != "#", l)
    vals = map(second, l)
    f.close()
    return vals

imrhovals = readsecond(fnim)
rerhovals = readsecond(fnre)
#os.unlink(fnim)
#os.unlink(fnre)

imAvals = readsecond(Afn)
reAvals = readsecond(reAfn)
os.unlink(reAfn)

imFvals = readsecond(Ffn)
reFvals = readsecond(reFfn)
os.unlink(reFfn)

# Calculate the selfenergy and the spectral function
fsigmaim = open('sigmaim.dat', 'w')
fsigmare = open('sigmare.dat', 'w')
f = open('spec.dat', 'w')

for i in range(len(xvals)):
  omega = xvals[i]
  # Raw Green's function
  Gold = 1j*imAvals[i] + reAvals[i]
  # F function
  F = 1j*imFvals[i] + reFvals[i]
  # Self-energy
  Sigma = J * F / Gold

  str = '%18.10f %18.15f\n' % (xvals[i], Sigma.imag)
  fsigmaim.write(str)
  str = '%18.10f %18.15f\n' % (xvals[i], Sigma.real)
  fsigmare.write(str)
  
  # G^0 (rho)
  G0 = -pi * (1j*imrhovals[i] + rerhovals[i])
  
  # Green's function
  G = 1.0/(1.0/G0 - Sigma)

  # Spectral density
  A = -1.0/pi * G.imag
  str = '%18.10f %18.15f\n' % (xvals[i], A)
  f.write(str)
  
fsigmaim.close()
fsigmare.close()
f.close()
