from sage.quadratic_forms.special_values import quadratic_L_function__exact as lfun
from sage.quadratic_forms.binary_qf import BinaryQF_reduced_representatives



#### BIG TODOS

## 1) get rid of redundant functions (ie functions already in sage)
## 2) mesh code with sage's BQF class. allow BQF related functions to take a BQF class as an argument as well
## 3) neaten thigns up., comment code, rename functions, add error handing. maybe some helpful print statements
## 4) makeCoeff cannot calculate for (0,0,0) fix this
#########################################################

#### NOTES

## 1) right now i have a couple ways of running this code
## one way is to run sage in terminal and import this file
# and execute things manually
# another is to do $ sage < file
# which is good step towards efficient troubleshooting since i dont
# have to manually enter commands every time
#  however when there are errors i cannot see
# print statements that i set up for troublshooting
# so instead i now run sage and use the load command
# when i use load("file.py") i get decimal remainders
# on coefficient values, which is puzzling
# i can fix this by doing load("file.sage") instead
# this leads me to suspect that my code has some
# inconsistencies between python and sage. im not sure
# what is going on. it seems i can press forward
# by using the last method i mentioned, but this
# anomaly should be kept in mind.


###
def mresidues(n, isprime = False):
    if isprime == True:
        A = list()
        for r in range(1,n):
            A.append(r)
        return A
    else:
        A = list()
        for r in range(1,n):
            if gcd(n, r) == 1:
                A.append(r)
        return A

#fundamental discriminant finder that Jake came up with, much easier to comprehend.

def fundDisc(N):
    sn = sign(N)
    N = abs(N)
    #strip out squares
    for x in range(2, int(floor(N ** 0.5)) + 1):
        if 0 == (N % (x**2)):
            return fundDisc(sn * N/(x**2))
    #re-sign the integer
    temp = sn * N
    if 1 == temp % 4:
        return temp
    return 4 * temp

def isFundDisc(N):
    if N == fundDisc(N):
        return True
    else:
        return False


#########################################################

##cohen's function, according to raum this can be calculated by sage... investigate this further.
def H(k1,N):
    #here we assign variables in a way that meshes with notation in McCarthy and Raum
    k = k1 + 1
    N = -N
    s = 2 - k
    #when N = 0 we calculate the value with the zeta func
    if N == 0:
        x = zeta((2 * s) - 1)
        return x
    #we check if N is 0 or 1 mod 4
    if N % 4 == 0 or N % 4 == 1:
        #find our kronecker character
        D0 = fundDisc(N)
        #find f
        f = isqrt(int(N / D0))
        #initialize our value
        hSum = 0
        #iterate through the divisors of f
        for d in range(1, f + 1):
            if f % d == 0:
                #when d divides f we add the appropriate term to the sum
                MOB = moebius(d)
                KRON = kronecker(D0, d)
                POW = d ** -s
                SIG = sigma(int(f / d),1 - (2 * s))
                #print('    adding ' + str(MOB) + ' * ' + str(KRON) + ' * ' + str(POW) + ' * ' + str(SIG) + ' to H func sum.')
                hSum += MOB * KRON * POW * SIG
        LF = lfun(s, D0)
        fCalc = hSum * LF
        return fCalc
    #otherwise the value is 0
    else:
        return 0


#########################################################

#this is where we throw everything together
def makeCoeff(B,k):
    if (B[0],B[1],B[2]) == (0,0,0):
        return 1
    if B.discriminant() > 0:
        B = B.reduced_form()
    #print(B)
    a=B[0]
    b=B[1]
    c=B[2]

    #calculate the discriminant
    disc = 4 * a * c - b * b
    
    #calculate the zeta functions
    z1 = zeta(1 - k)
    z2 = zeta(3 - (2 * k))

    #make a list of nonzero coefficients of our BQF
    list = []
    for i in [a, b, c]:
        if i != 0:
            list.append(i)
            
    sm = 0
    
    #iterate through divisors of our BQF (this forms the core sum)
    for d in range(1, abs(min(list)) + 1):
        if a % d == 0 and b % d == 0 and c % d == 0:
            #add term as given in McCarthy for each divisor d
            POW = (d ** (k-1))
            HFUN = H(k-1, disc / (d ** 2))
            sm += POW * HFUN
    #finish off the calculation by multiplying by 2/(Z(1 - k) * Z(3 - 32))        
    value = 2 * sm / (z1 * z2) 
    return value

#tada!


    
#print('The coefficient indexed by (' + str(a) + ', ' + str(b) + ', ' + str(c) +
#') for the Siegel Eisenstein series of weight ' + str(k) +   ' is ' + str(makeCoeff(a,b,c,k)))


##########################################



def giveReps(D):
    A = BinaryQF_reduced_representatives(D)
    return A


#iterate through some discriminants and their reduced forms to make some coeffs
#todo : change index from bqf to tuple, fix coefficients having a floating ".0"
def genJson(k, N):
    A = OrderedDict()
    for D in range(N):
        if isFundDisc(-D) == True:
            B = dict()
            for Y in giveReps(-D):
                B[str(y)] = str(makeCoeff(Y,k))
            A[str(-D)] = B
    return A

#########################

def BQFsplitter(B):
    a=B[0]
    b=B[1]
    c=B[2]
    D = list()
    for n in range(0,a+1):
        a1 = a - n
        for m in range(0,c+1):
            c1 = c - m
            ran = floor(2*sqrt(a1 * c1))
            for k in range(-int(ran),int(ran)+1):
                S1cand = BinaryQF(a1,k,c1)
                b2 = b-k
                S2cand = BinaryQF(n,b2,m)
                if (b2)**2 <= 4*m*n:
                    D.append((S1cand,S2cand))
    return D


#########################

def bracket(A,S):
    return S.transpose()*A*S

def matToBQF(A):
    if A[0,1] != A[1,0]:
        pass
        #raise error 
    a = int(A[0,0])
    b = int(2*A[0,1])
    c = int(A[1,1])
    return BinaryQF(a,b,c)

#############
#Class stuff#
#############


### Class stuff not really good for multiplying forms
#thinking of making a new structure where we name a form a 
# 1) "base" or "eisen" type in which case coeffs are gotten from the eisenceoffcalc
# 2) "add" type where we add the coeffs from the two forms that make it up
# 3) "mult" type where we get the coeffs by using a multiply routine on the two constituent forms. (using BQFsplitter)
# 4) "scale" type where we get the coeffs by a simple scalar multiplication

class SPMF: 
    def __init__(self, k):
        self.weight = k
        self.type = 'eisen'
        self.constituents = (self,1)

        
    def coeff(self,B):
        if self.type == 'eisen':
            return makeCoeff(B, self.weight)
        if self.type == 'mult':
            coeff = 0
            F1 = self.constituents[0]
            F2 = self.constituents[1]
            for Q in BQFsplitter(B):
                Q1 = Q[0]
                Q2 = Q[1]
                coeff += (F1.coeff(Q1) * F2.coeff(Q2))
            return coeff    
        if self.type == 'add':
            F1 = self.constituents[0]
            F2 = self.constituents[1]
            coeff = F1.coeff(B) + F2.coeff(B)
            return coeff
        if self.type == 'scale':
            F1 = self.constituents[0]
            c1 = self.constituents[1]
            coeff = c1 * F1.coeff(B)
            return coeff
        if self.type == 'twist':
            F1 = self.constituents[0]
            p = self.constituents[1]
            a = B[0]
            b = B[1]
            c = B[2]
            #case (i)
            if ((not b % p) == 0 and a % p**4 == 0):
                P_1 = p**(1 - F1.weight) * kronecker(b,p)
                P_2 = 0
                A=matrix([[a,Rational(b/2)],[Rational(b/2),c]])
                for r in mresidues(p, isprime = True):
                    S = matrix([[1,Rational(-r/p)],[0,p]])
                    #print('mark1')
                    newS = matToBQF(bracket(A,S))
                    print(newS)
                    P_2 += kronecker(r,p) * F1.coeff(newS)
                    #print('mark3')
                P = P_1 * P_2
                return P
            else:
                return 'conditionsnotmet'

def multiplyForms(E,F):
    A = SPMF(E.weight + F.weight)
    A.type = 'mult'
    A.constituents = (E,F)
    return A

def addForms(E,F):
    A = SPMF(E.weight)
    A.type = 'add'
    A.constituents = (E,F)
    return A

def scaleForms(E,c):
    A = SPMF(E.weight)
    A.type = 'scale'
    A.constituents = (E,c)
    return A

def twistForms(E,p):
    A = SPMF(E.weight)
    A.type = 'twist'
    A.constituents = (E,p)
    return A