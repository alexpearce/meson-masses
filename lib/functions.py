from math import fabs

""" Returns -1, 0, +1 for x < =, x == 0 and x > 0 respectively """
sgn = lambda x : cmp(x,0)

def sign(a, b):
  """Returns |a| if b >= 0, -|a| otherwise"""
  return fabs(a) if b >= 0 else -fabs(a)