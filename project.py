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
from math import pow, factorial, fabs

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
  x_3  = (x*x*x) / 3.0
  
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

def c_k_arr(k):
  """
  Computes the first k+2 c_k coefficients for integer k > 1
  Interestingly, this is correct for c_k k > 6, whereas c_k(k) is not
  Good for about 14 dp up to c_k k > 9, not tested above this
  """
  c = range(k + 2)
  c[0] = 1.0
  c[1] = 3.0*5.0 / 216.0
  for i in range(2, k+2):
    six_i = 6.0*i
    numerator   = (six_i - 5)*(six_i - 3)*(six_i - 1)
    denominator = 216.0 * i * ((2.0*i) - 1)
    
    c[i] = c[i-1] * (numerator/denominator)
  
  return c

def airy_two(x, n = 10):
  """Computes the second approximation of Ai(x)"""
  zeta      = (2.0/3.0) * pow(x, (3.0/2.0))
  prefactor = 0.5 * pow(pi, -0.5) * pow(x, -0.25) * exp(-zeta)
  
  s = zeros(n)
  for i in range(n):
    s[i] = pow(-1.0, i) * c_k(int(i)) * pow(zeta, -i)
  
  return prefactor * sum(s)
  
def airy_three(x, n = 3):
  """Computes the third approximation of Ai(x)"""
  mod_x = fabs(x)
  zeta  = (2.0/3.0) * pow(mod_x, (1.5))
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
  
def find_roots(f, a, b, h = 0.01):
  """
  Returns tuples of (x_1, x_2), where between x_1 and x_2 Ai(x) = 0,
  for Ai(x) between a and b, 'walking' along the function in steps of h
  """
  rng = arange(a, b, h)
  # All the zeros are for x < 0 so use the third approximation for now
  ai_x = [f(i) for i in rng]
  
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

### Roots ###

def roots():
  """
  All the roots between -2 and -15
  """
  # Use the correct approximation in each range
  a_one_roots = find_roots(airy_one, -2, -7, -0.01)
  # Join the two approximations correctly, and not near a root
  # (by inspection)
  a_two_roots = find_roots(airy_three, -7, -15.01, -0.01)
  
  rts = []
  
  for root in a_one_roots:
    rts.append(ridders_method(airy_one, root[0], root[1], 1e-10))
  for root in a_two_roots:
    rts.append(ridders_method(airy_three, root[0], root[1], 1e-10))
    
  return rts

def plot_points():
  """
  Plots a pretty graph of the numerical roots of Airy's function
  against the experimental meson masses. Also finds and prints the
  slope and intersection of the LOBF for both sets of points.
  """
  # Take the minus sign as the roots are negative (and it says we should in E = -b/a)
  energies = [-rt for rt in roots()]

  b_masses = [9.46, 10.02, 10.35, 10.57, 10.86, 11.02]
  c_masses = [3.10, 3.69, 4.04]

  b_energies = energies[:6]
  c_energies = energies[:3]

  (b_slope, b_intercept) = polyfit(b_energies, b_masses, 1)
  (c_slope, c_intercept) = polyfit(c_energies, c_masses, 1)

  # Theory bottom bass is 4.19 GeV
  print "b mass: {}, b slope: {}".format(b_intercept/2, b_slope)
  # Theory charm mass is 1.29 GeV
  print "c mass: {}, c slope: {}".format(c_intercept/2, c_slope)

  # Plot the numerical data vs the experimental
  plb.scatter(b_energies, b_masses)
  plb.scatter(c_energies, c_masses)

  # Plot the fitted curves
  rng = arange(0, 11, 1)
  plb.plot(rng, [polyval([b_slope, b_intercept], x) for x in rng])
  plb.plot(rng, [polyval([c_slope, c_intercept], x) for x in rng])

  # Set the y-axi range and display a grid
  plb.yticks(arange(0, 12))
  plb.xticks(arange(0, 10))
  plb.ylim([0, 12])
  plb.xlim([0, 10])
  plb.grid(True)

  # Axes labels
  plb.xlabel("Numerical")
  plb.ylabel("Experimental")
  plb.title("s-state meson masses")

  plb.show()