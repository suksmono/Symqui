from sympy import *
import itertools as itr
import numpy as np


#from . import symq
'''
----------
* symqui *
----------
A library of symbolic computing for QUBO problems
qubo: quadrature binary optimization

----------
howto use:
----------
# import library
    >>import symque as sq
#define object, eg.               
    >> qb=sq.symqx()
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
def eq_standardized(Hk,sqChar):
    '''
    change all input variables into a standard variables
    eg. E = x1 + x2*x4  -> E=s1 + s2*s3 
    assume EQU is already symbolic
    '''
    # get a list of all variables in the expression
    # get the symbols ordered  [x0, x1, x2, ...]
    symbSet=tuple(ordered(Hk.free_symbols))
    symOrg=[]
    symStd=[]
    idx=0
    for tx in symbSet:
        symOrg.append(tx)
        if (sqChar.upper()=='S'):
            tsq=symbols('s%d'%idx)
        else:
            tsq=symbols('q%d'%idx)
        symStd.append(tsq)
        # Hamiltonian with standardized variables
        Hk=Hk.subs({tx:tsq})
        '''
        symOrg[x0, x1, ..., xk]
        symStd[s0, s1, ..., sk]
        '''
        #print('idx substitution->',idx)
        idx=idx+1
        #
    return symOrg, symStd, Hk
##
def simplify_sq_squares(Hkb, tSq):
    '''
    simplify (binary symbolic) function by assigning all
    squared variables with one: qi**2 <- qi, si**2 <-1  
    where bi is a binary variable, either qi={0,1} or si{-1,+1}
    '''
    print('Simplifying bi**2 <-1 ')
    # get all free symbols in Hkb
    bb=Hkb.free_symbols
    # substitute
    for tb in bb:
        if(tSq.upper()=='S'):
            Hkb= Hkb.subs( {tb**2:1}) 
            Hkb= Hkb.subs( {tb**3:tb}) 
            #print('Processing qi^2->1:',ts)
        else:
            Hkb= Hkb.subs( {tb**2:tb}) 
            Hkb= Hkb.subs( {tb**3:tb}) 
         
    return(simplify(Hkb))
##
def H2sub(x1,x2,y,d12):
    '''
    input:  x1, x2 of x1*x2 product
            y is var result x1*x2 <-- y
            d12 compensation factor
    output:two-body polynomial compensation term
    '''
    return(d12*(3*y+x1*x2-2*x1*y-2*x2*y))

#
def boole_reduce(Hkq,delta):
    '''
    input:  k-body hamiltonian in q-domain Hkq
    output: 2-body hamiltonian in q-domain H2s
    stages: Hkq -> H2q
    '''
    # list all involved variables
    # define higher order set >2 
    hoT=set() 
    # collection all variables in Hk
    toT=set()
    # convert Hks->Hkq
    #Hkq=Hks2Hkq(Hks)
    tSet= set(Hkq.args)
    for tt in tSet:
        # get all variables in a term
        ta=tt.free_symbols # obtain var only
        toT=toT.union(ta)
        #print('ta->',ta)
        if( len(ta)>2 ):
            hoT=hoT.union(ta)
    # knowing vadiables in high order term, now construct substitution list
    # the index of ancillary variable should start after number of var in Hk
    #nvTO=len(toT) # number of total variables in Hk
    nvHO=len(hoT) # number of variables involved higher order terms in Hk 
    #
    print('hoT variables->',hoT)
    qPair=gen_subtitution_pairs(hoT,nvHO)
    # ##debug##
    #print('*all substitution pairs',qPair)
    #delta=100
    d=delta; #2*vMax 
    '''
    # do substitution iteratively 
    '''
    H2q=Hkq
    for tx in qPair:
        if (is_ho_equation(H2q)==True):   
            #print(tx)
            # do substitition
            ## do substitution only to high-order terms
            print(tx[0],'*' , tx[1], '->',tx[2] )
            H2q = H2q.subs({tx[0]*tx[1]:tx[2]} ) \
                    + H2sub(tx[0], tx[1],tx[2],d)
    
            H2q=simplify(H2q)
            print('H2q->', H2q)
            print('\n')
    #
    return expand(H2q) 
##
##
def gen_subtitution_pairs(hoT,idxStart):
    '''
    input: hoT-a set of high order variables 
           idxStart -start index of ancillary variables
    output: list of substitution pair, [q1,q2,q3] means:  (q1*q2)<--q3
    '''
    #
    lsubP=[]
    #print('hoT->',hoT)
    hoTlist=list(ordered(hoT))
    #print('hoT-list ->',hoTlist)
    #cp=set(itr.combinations(hoTlist,2))
    cp=list(itr.combinations(hoTlist,2))
    print('cp->',cp)
    idx=idxStart
    for tz in cp:
        ttz=list(tz)
        # add ancillary variable to var pair->triple
        ttz.append(symbols('q%d'%idx))
        idx=idx+1
        lsubP.append(ttz)
    ##
    #print('*all substitution pairs',lsubP)
    return lsubP
###
####################################################
def is_ho_equation(tEq):
    '''
    check, whether an Equation is of high order, i.e., 
    it has a term of order > 2
    input : symbolic Equation, eg. 1+q0*q1, 2+s3**2+s0*s1*s2 ...
    output: True/False
    '''
    oTerms=tEq.as_ordered_terms()
    #print('term', oTerms,'->free symbols:', oTerms[0].free_symbols)
    if len( oTerms[0].free_symbols )>2:    #order: decreasing
        return True
    else:
        return False
####################################################  
#
############################################################

def get_qs_coeffs(H2qs,qsType):
    '''
    get Ising or qubo coefficients
    input:  > H2qs  :a 2-body hamiltonian
            > qsType: problem type=> 'q'/0: qubo, 's'/1: ising
    output: qs_coeffs {hi, Jij}
    assume the symbols are
    s0, s1, ...., s_(NQ-1) or q0, q1, ...., q_(NQ-1) 
    '''
    NQ=len(H2qs.free_symbols)
    hi=np.zeros(NQ)
    Jij=np.zeros((NQ,NQ))
    dc=H2qs.as_coefficients_dict()
    # list all symbols: s0, s1, ...
    if(qsType.upper()=='S'):
        qs=symbols('s0:%d'%NQ)
    else:
        qs=symbols('q0:%d'%NQ)

    # extract b
    b=dc[1]
    # extract hi
    for m in range(NQ):
        hi[m]=dc[qs[m]];
    # extract Jij
    for m in range(NQ):
        for n in range(m+1,NQ):
            print('m=',m,'n=',n)
            Jij[m,n]=dc[qs[m]*qs[n]]
    # return the results
    return(b, hi, Jij)
############################################################
class symqx:
    '''
    qubo-cefficient extractor
    '''
    
    def __init__(self):
        self.eq_str=''
        self.Ekx=''
        # list of original variables
        self.var_original=[]
        # list of standardized variables
        self.var_standard=[]

        self.Ekq=''
        #self.Hks=''
        self.Hkq=''
        #self.H2s=''
        self.H2q=''
        self.verbose=True
        self.delta=100
        self.qubo_b=0
        self.qubo_hi=[]
        self.qubo_Jij=[]
        self.wildqat_Jij=[] # wildqat style Jij, embed hi in Jij 
        self.wildqat_flag=True
    def get_qubo_coeffs(self):
        '''
        a problem given in string format is stored in eq_str
        '''
        self.Ekx=sympify(self.eq_str,evaluate=False)
        print('Ekx=', self.Ekx)
        self.var_original,self.var_standard,self.Ekq = eq_standardized(self.Ekx,'q')
        print('Ekq=', self.Ekq)

        # simplifying qi**2 <- 1
        self.Hkq=expand(simplify_sq_squares(expand(self.Ekq),'q'))
        if self.verbose:            
            print('Simplified Hkq=', self.Hkq)
        # Ekq-> E2q
        self.H2q =boole_reduce(self.Hkq,self.delta)
        #self.H2q=H2q
        if self.verbose:            
            print('H2q=', self.H2q)
        #return 0
        self.qubo_b, self.qubo_hi, self.qubo_Jij=get_qs_coeffs(self.H2q,'q')
        if self.verbose:            
            print('Jij=', self.qubo_Jij)
        #
        self.wildqat_Jij=self.qubo_Jij
        for m in range(0,len(self.qubo_hi)):
            self.wildqat_Jij[m][m]=self.qubo_hi[m]
        if self.verbose:            
            print('Jij=', self.wildqat_Jij)
        #
        if self.wildqat_flag:            
            return self.wildqat_Jij.tolist()
        
        # 