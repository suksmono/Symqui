'''
----------
* symque *
----------
A library of symbolic computing for qubo problems
qubo: quadrature binary optimization

----------
howto use:
----------
# import library
    >>
#define object, eg.               
    >> qb = symqx()
#enter the problem in string, eg. 
    >> qb.eq_str = '2*(1+q0*q1*q2*q3)'  
#extract qubo coefficients
    >> J=qb.get_qubo_coeffs()
#define solver
    >>import wildqat as wq
    >>a = wq.opt()
#enter qubo coeffs into solver
    >>a.qubo=J
#solve the problem
    >>a.solve()
'''
#
from sympy import *
#import itertools as itr
#import numpy as np
#import neal
#

class symqx:
    '''
    qubo-cefficient extractor
    '''
    
    def __init__(self):
        self.eq_str=''
        self.Ekx=''
        #self.Ekq=''
        #self.Hks=''
        #self.Hkq=''
        #self.H2s=''
        #self.H2q=''
    def get_qubo_coeffs(self):
        '''
        a problem given in string format is stored in eq_str
        '''
        self.Ekx=sympify(self.eq_str,evaluate=False)
        print('Ekx=', self.Ekx)

        