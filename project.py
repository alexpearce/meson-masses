#!/usr/local/bin/python

# Import 3rd party modules
import pylab as plb
from numpy import *
# Import Airy's functions
# http://docs.scipy.org/doc/scipy/reference/generated/scipy.special.airy.html
from scipy.special import airy

# Import our modules
# from lib.module_name import *

def f(x, n = 100):
  """Returns f(x) summed to n terms"""
  f = zeros(n)
  
  f[0] = 1
  
  x_3 = (x*x*x) / 3.0
  for i in range(1, n): 
    f[i] = (f[i-1]*x_3) / ((3.0*i*i) - i)
  
  return sum(f)  

def g(x, n = 100):
  """Returns g(x) summed to n terms"""
  g = zeros(n)
  
  g[0] = x
  x_3 = (x*x*x) / 3.0
  
  for i in range(1, n):
    g[i] = (g[i-1]*x_3) / ((3.0*i*i) + i)
    
  return sum(g)

def airy_one(x, n = 100):
  """
  Computes the first approximation of Ai(x)
  using Ai(x) = a*f(x) - b*g(x)
  """
  a = 0.3550280538
  b = 0.2588194037
  
  return (a*f(x, n)) - (b*g(x, n))
  
def compare_airy_one(a, b, n = 100):
  """
  Compares the first approximation of Ai(x)
  with numpy's Ai(x), in the range [a, b), b > a,
  using n terms in the approximation sum.
  """
  
  # Step size of 0.01
  x = arange(a, b, 0.01)
  
  theory = [airy(y)[0] for y in x]
  # Compute n terms in the sum
  experiment = [airy_one(y, n) for y in x]

  plb.plot(x, theory)
  plb.plot(x, experiment)
  plb.show()