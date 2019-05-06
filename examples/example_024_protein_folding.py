# -*- coding: utf-8 -*-
"""
Created on Fri May  3 12:56:26 2019
@author: Suksmono
Problem: protein folding
"""
'''
# setup symqui
'''
import symqui as sq
qb=sq.symqx()
qb.delta=10
'''
enter the equation to solve
'''
qb.eq_str = '-1-4*q3+9*q1*q3+9*q2*q3-16*q1*q2*q3'

'''
# solve using wildqat
'''
import wildqat as wq
a = wq.opt()
#enter qubo coeffs into solver
a.qubo= qb.get_qubo_coeffs() 
#solve the problem
print('\nSolution:',a.sa())
