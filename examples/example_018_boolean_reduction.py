# -*- coding: utf-8 -*-
"""
Created on Fri May  3 12:56:26 2019
@author: Suksmono
Problem: boolean reduction
"""

'''
#-----------------------------------
# setup symqui
#-----------------------------------
'''
import symqui as sq
qb=sq.symqx()
qb.delta=0.5
qb.eq_str = 'q0-q0*q1*q2' #enter the function as string

'''
#-----------------------------------
# solve using wildqat
#-----------------------------------
'''
import wildqat as wq
a = wq.opt()
#enter qubo coeffs into solver
a.qubo= qb.get_qubo_coeffs() 
#solve the problem
print('\nSolution:',a.sa())
