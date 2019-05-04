# -*- coding: utf-8 -*-
"""
Created on Fri May  3 12:56:26 2019
@author: Suksmono
Problem: prime factorization
"""
'''
# setup symqui
'''
import symqui as sq
qb=sq.symqx()
qb.delta=128
'''
enter the equation to solve
'''
qb.eq_str = '(15-(1+2*x0+4*x1)*(1+2*x2))**2'

'''
# solve using wildqat
'''
import wildqat as wq
a = wq.opt()
#enter qubo coeffs into solver
a.qubo= qb.get_qubo_coeffs() 
#solve the problem
print(a.sa())
