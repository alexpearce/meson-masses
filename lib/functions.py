from math import fabs, sqrt

""" Returns -1, 0, +1 for x < =, x == 0 and x > 0 respectively """
sgn = lambda x : cmp(x,0)

def sign(a, b):
  """Returns |a| if b >= 0, -|a| otherwise"""
  return fabs(a) if b >= 0 else -fabs(a)

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