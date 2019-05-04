# examples-01
# k to 2 body

import symque as sq

# create an object qsymaqc 
mypb=sq.symq
# fill EQ with expression in STRING format
mypb.EQ='x0 - x0*x1*x2'

# extract qubo/ising coefficients
hi, Jij = mypb.extract_qs_coeffs_from_string() 
J=Jij.tolist()
for m in range(0,len(hi)):
    J[m][m]=hi[m]
## to be processed by neal, wildqat, qbsolve, dwave, or other quantum machines
#print('\nhi=',mypb.ising_hi)
#print('\nJij=',mypb.ising_Jij)

###
'''
Solve the problem using wildqat
'''
print('\nSolve using wildqat ... !')
import wildqat as wq
b = wq.opt()
# fill in qubo coeffs
b.qubo = J#Jij.tolist() 
print(b.sa())