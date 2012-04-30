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
from lib.functions import *



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
  s_1, s_2 = zeros(n), zeros(n)
  for i in range(n):
    s_1[i] = pow(-1.0, i) * c_k(2*i) * pow(zeta, -2.0*i)
    s_2[i] = pow(-1.0, i) * c_k((2*i) + 1) * pow(zeta, -(2*i) - 1)
  
  return prefactor * (sin(trig_arg)*sum(s_1) - cos(trig_arg)*sum(s_2))

def nice_example():
  """A nice example of how the three approximations work together"""
  # The ranges across which each approximation will act
  a_one_rng   = linspace(-15, 11, 1000)
  a_two_rng   = linspace(1,   20, 1000)
  a_three_rng = linspace(-20, -1, 1000)

  a_one   = [airy_one(x)   for x in a_one_rng]
  a_two   = [airy_two(x)   for x in a_two_rng]
  a_three = [airy_three(x) for x in a_three_rng]

  plb.plot(a_one_rng, a_one)
  plb.plot(a_two_rng, a_two)
  plb.plot(a_three_rng, a_three)
  plb.show()
  
def find_roots(a, b, h = 0.01):
  """
  Returns tuples of (x_1, x_2), where between x_1 and x_2 Ai(x) = 0,
  for Ai(x) between a and b, 'walking' along the function in steps of h
  """
  rng = arange(a, b, h)
  # All the zeros are for x < 0 so use the third approximation for now
  ai_x = [airy_one(i) for i in rng]
  
  prev_sign = 0
  
  sign_changes = []

  for idx, val in enumerate(ai_x):
    # 'Walk' along the function watching for sign changes
    new_sign = sgn(val)
    
    if (prev_sign + new_sign) == 0:
      # -1 + 1 = 0, so we've got a sign change
      # tuple of the current x value and the previous one
      sign_changes.append((rng[idx], rng[idx - 1]))
    
    prev_sign = new_sign
    
  return sign_changes

def ridders_method(f, x1, x2, xacc = 0.001):
  """
  Finds a root of f between x1 and x2 to precision xacc
  Translated to Python from Numerical Recipes 2007, p453
  http://apps.nrbook.com/empanel/index.html#pg=453
  """
  iterations = 50
  f1 = f(x1)
  f2 = f(x2)
  
  if (sgn(f1) + sgn(f2)) == 0:
    # A highlight unlikely root value
    ans = -9e15
    for i in range(0, iterations):
      # Midpoint
      x3 = 0.5*(x1 + x2)
      f3 = f(x3)
      denominator = sqrt((f3*f3) - (f1*f2))
      # Check for 1/0
      if denominator == 0: return ans
      x4 = x3 + ((x3 - x1)*(sgn(f1 - f2)*f3)/denominator)
      if fabs(x4 - ans) <= xacc: return ans
      ans = x4
      f4 = f(ans)
      if f4 == 0: return ans
      
      if sign(f3, f4) != f3:
        # Root is in (x3, x4)
        x1 = x3
        f1 = f3
        x2 = ans
        f2 = f4
      elif sign(f1, f4) != f1:
        # Root is in (x1, x4)
        x2 = ans
        f2 = f4
      elif sign(f2, f4) != f2:
        # Root is in (x2, x4)
        x1 = ans
        f1 = f4
      else:
        exit("Something funny's going on.")
      
      if fabs(x2 - x1) <= xacc: return ans
    
    exit("Ridders' could not acheive desired accuracy.")
  else:
    if f1 == 0: return x1
    if f2 == 0: return x2
    exit("No root found between x1 and x2.")

rts = find_roots(-8, -3)
for tup in rts:
  print ridders_method(airy_one, tup[1], tup[0], 1e-7)