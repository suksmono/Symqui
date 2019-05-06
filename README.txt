#-------------------------------------------------------
# SYMQUI: (SYMbolic Quantum User Interface)
#-------------------------------------------------------
# A symbolic quantum computing interface to qubo-solver
# given a string expression of binary {0,1} variable,
# symque will present the qubo-coefficient to be solved
# by a qubo solver, such as wildqat or qbsolve 
#
# Work flow:
# ---------
#
# Eq_string -> Ekx -> Ekq -> E2q -> {hi, Jij}
#
# {string_equation}->|QSymbolic|->{qubo_coeffs}->|qubo_solver|->{Results}
#
# legend: {... data ...}, | ... process ...|   


#---------------
# how to install
#---------------

    >>pip install symqui

#---------------
# how to use:
#---------------
# import the library
    >>import symqui as sq
# define a qubo extractor object, eg.               
    >>qb=sq.symqx()
# enter the problem in string, eg. 
    >>qb.eq_str = '2*(1+q0*q1*q2*q3)'  
# extract qubo coefficients
    >>J=qb.get_qubo_coeffs()
# define solver
    >>import wildqat as wq
    >>a = wq.opt()
# enter qubo coeffs into solver
    >>a.qubo=J
# solve the problem
    >>a.solve()
                                                                                                                                                               