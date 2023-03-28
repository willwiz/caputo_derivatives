# File: model.pyx
'''
This file implements the Holzapfel Gasser Ogden model 2013 version, updated in 2015 in
    Holzapfel et al. Modelling non-symmetric collagen ... J. Royal Soc. Interface. 2015

Jan 15, 2019
By: Will Zhang
'''

'''
HOG-2013 Model
Parameter List
mu      :   Matrix modulus
k1      :   Collagen modulus
k2      :   Collagen exponent
e1      :   Elastin modulus
e2      :   Elastin exponent
s1      :   Smooth muscle exponent
s2      :   Smooth muscle exponent

Vector
alpha   :   Fiber split
kop     :   splay parameters
kip     :   splay parameters

Inputs:
C11
C12
C21
C22

Outputs:
S11
S12
S21
S22

'''


from numpy import exp, ndarray, array, cos, sin


identity2D = array([[1,0,0,1]],dtype=float)

class hog20132D:
  def __init__(self, pars: ndarray, alpha:float, kip:float, kop:float) -> None:
    # unpack parameters
    self.mu = pars[0]
    self.k1 = pars[1]
    self.k2 = pars[2]
    self.e1 = pars[3]
    self.e2 = pars[4]
    self.s1 = pars[5]
    self.s2 = pars[6]
    # Compute fiber info
    self.A = 2.0*kop*kip
    self.B = 2.0*kop*(1.0-2.0*kip)
    self.C = 1.0 - 3.0*self.A - self.B
    # Compute the othogonal fiber vector
    cosa = cos(alpha)
    sina = sin(alpha)
    self.m4 = array([cosa*cosa,  cosa*sina,  cosa*sina, sina*sina])
    self.m6 = array([cosa*cosa, -cosa*sina, -cosa*sina, sina*sina])
    self.H4 = self.A*identity2D + self.B*self.m4
    self.H6 = self.A*identity2D + self.B*self.m6
    # print(self.H4 + self.H6)
    self.He = array([1,0,0,0], dtype=float)
    self.Hs = array([0,0,0,1], dtype=float)
  def stress(self, args : ndarray) -> ndarray:
    det = (args[0]*args[3] - args[1]*args[1])
    I_n = 1.0 / det
    Cinv = array([args[3], -args[1], -args[2], args[0]]) * I_n
    I_1 = args[0] + args[3] + I_n
    I_4 = self.m4 @ args
    I_6 = self.m6 @ args
    E6  = self.A*I_1 + self.C*I_n - 1.0
    E4  = E6 + self.B*I_4
    E6  = E6 + self.B*I_6
    dWd4 = self.k1*E4*exp(self.k2*E4*E4)
    dWd6 = self.k1*E6*exp(self.k2*E6*E6)
    # print(I_4, I_6, dWd4, dWd6)
    Ee = args[0] - 1.0
    Es = args[3] - 1.0
    dWde = self.e1*Ee*exp(self.e2*Ee*Ee)
    dWds = self.s1*Es*exp(self.s2*Es*Es)
    p = (self.mu + self.C * dWd4 + self.C * dWd6) * I_n
    val = self.mu * identity2D + dWd4*self.H4 + dWd6*self.H6 + dWde*self.He + dWds*self.Hs
    return val - p*Cinv
  def stress_parts(self, args : ndarray) -> ndarray:
    det = (args[0]*args[3] - args[1]*args[1])
    I_n = 1.0 / det
    Cinv = array([args[3], -args[1], -args[2], args[0]]) * I_n
    I_1 = args[0] + args[3] + I_n
    I_4 = self.m4 @ args
    I_6 = self.m6 @ args
    E6  = self.A*I_1 + self.C*I_n - 1.0
    E4  = E6 + self.B*I_4
    E6  = E6 + self.B*I_6
    dWd4 = self.k1*E4*exp(self.k2*E4*E4)
    dWd6 = self.k1*E6*exp(self.k2*E6*E6)
    # print(I_4, I_6, dWd4, dWd6)
    Ee = args[0] - 1.0
    Es = args[3] - 1.0
    dWde = self.e1*Ee*exp(self.e2*Ee*Ee)
    dWds = self.s1*Es*exp(self.s2*Es*Es)
    # p = (self.mu + self.C * dWd4 + self.C * dWd6)
    # val = self.mu * identity2D + dWd4*self.H4 + dWd6*self.H6 + dWde*self.He + dWds*self.Hs
    return self.mu * identity2D, dWd4*self.H4 + dWd6*self.H6, dWde*self.He, dWds*self.Hs, array([self.mu]), array([self.C * dWd4 + self.C * dWd6]), I_n*Cinv









