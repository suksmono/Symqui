# -*- coding: utf-8 -*-
"""
Created on Fri May  3 12:56:26 2019
@author: Suksmono
Simplest example of symque
"""
'''
# setup symque
'''
import symqui as sq
qb=sq.symqx()
qb.eq_str = '2*(1+q0*q1*q2*q3)'

'''
# solve using wildqat
'''
import wildqat as wq
a = wq.opt()
#enter qubo coeffs into solver
a.qubo= qb.get_qubo_coeffs() 
#solve the problem
print('\nSolution:',a.sa())
