#!/usr/local/bin/python

# Import 3rd party modules
# Plotting
import pylab as plb
# zeros, arange etc.
from numpy import *
# Import Airy's functions
# http://docs.scipy.org/doc/scipy/reference/generated/scipy.special.airy.html
from scipy.special import airy
# Import a few maths functions
from math import pow, factorial

# Import our modules
# from lib.module_name import *



## The First Approximation ##

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
  It is satisfactorily stable for (approx.) -14 < x < 10
  """
  a = 0.3550280538
  b = 0.2588194037
  
  return (a*f(x, n)) - (b*g(x, n))
  
def compare_airy_one(a = -14, b = 10, n = 100):
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



## The Second Approximation ##

def c_k(k):
  """Computes the coefficient c_k for integer k > 0"""
  if isinstance(k, int) == False or k < 0: exit("k is not a positive integer")
  if k == 0: return 1.0
  
  two_k = 2*k
  limit = (3*two_k) - 1
  
  # Multiply together all the elements in `range`
  numerator   = multiply.reduce(range(two_k + 1, limit + 1, 2))
  denominator = pow(216, k) * factorial(k)
  
  return (numerator / denominator)

def airy_two(x, n = 10):
  """Computes the second approximation of Ai(x)"""
  zeta = (2.0/3.0) * pow(x, (3.0/2.0))
  prefactor = 0.5 * pow(pi, -0.5) * pow(x, -0.25) * exp(-zeta)
  
  s = zeros(n)
  for i in range(n):
    s[i] = pow(-1.0, i) * c_k(int(i)) * pow(zeta, -i)
  
  return prefactor * sum(s)
  
def airy_three(x, n = 2):
  """Computes the third approximation of Ai(x)"""
  mod_x = fabs(x)
  zeta = (2.0/3.0) * pow(mod_x, (1.5))
  prefactor = pow(pi, -0.5) * pow(mod_x, -0.25)
  
  trig_arg = zeta + (pi/4.0)
  
  # Compute the first and second sum
  s_1 = zeros(n)
  s_2 = zeros(n)
  for i in range(n):
    s_1[i] = pow(-1.0, i) * c_k(2*i) * pow(zeta, -2.0*i)
    s_2[i] = pow(-1.0, i) * c_k((2*i) + 1) * pow(zeta, -(2*i) - 1)
  
  return prefactor * (sin(trig_arg)*sum(s_1) - cos(trig_arg)*sum(s_2))

def nice_example():
  """A nice example of how the three approximations work together"""
  a_one_rng = linspace(-15.5, 11, 1000)
  a_two_rng = linspace(5, 20, 100)
  a_three_rng = linspace(-20, -5, 1000)
  a_one   = [airy_one(x) for x in a_one_rng]
  a_two   = [airy_two(x) for x in a_two_rng]
  a_three = [airy_three(x) for x in a_three_rng]

  plb.plot(a_one_rng, a_one)
  plb.plot(a_two_rng, a_two)
  plb.plot(a_three_rng, a_three)
  plb.show()