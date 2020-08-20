from sage.functions.transcendental import zeta
from sage.quadratic_forms.special_values import quadratic_L_function__exact as lfun
from sage.quadratic_forms.binary_qf import BinaryQF_reduced_representatives
from sage.arith.misc import jacobi_symbol
import numpy as np
import math
from fractions import Fraction
from collections import OrderedDict


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

#here we have a simple prime checker
#utilzes theorem that says n composite ==> n has prime divisor less than sqrt(n)
#so the contrapositive says if n has no prime divisor less than sqrt(n) then n is prime

def isPrime(n):
    if (n < 2) :
        return False
    for i in range(2, n + 1) :
        if (i * i <= n and n % i == 0) :
            return False
    return True

def leg(a,b):
    return jacobi_symbol(a,b)


#########################################################

#here we sum over the divisors of n raised to the kth power
#this is known as the divisor function

def sigmaDivisor(k, n):
    y = 0
    #could cut this to only go to n+1/2 if n gets big
    #~~it probably will~~
    for x in range(1,int(n+1)):
        #we check if x is a divisor of n
        if n % x == 0:
            #if it is a divisor we add its kth power to the value
            y += x ** k
    return int(y)


#########################################################

#found this function online, but it seems to work and the code makes sense so we use it

def mobius(N) :
    # Base Case
    if (N == 1) :
        return 1
    # For a prime factor i
    # check if i^2 is also
    # a factor.
    p = 0
    for i in range(1, N + 1) :
        if (N % i == 0 and
                isPrime(i)) :
            # Check if N is
            # divisible by i^2
            if (N % (i * i) == 0) :
                return 0
            else :
                # i occurs only once,
                # increase f
                p = p + 1
    # All prime factors are
    # contained only once
    # Return 1 if p is even
    # else -1
    if(p % 2 != 0) :
        return -1
    else :
        return 1


#########################################################

def gcd(a, b):
    while(b):
        a, b = b, a%b
    return a

def divides(a,b):
    if b % a == 0:
        return True
    else:
        return False

def maxdivides(p,k,b):
    if b % p**k == 0 and b % p**(k+1) != 0:
        return True
    else:
        return False

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

#########################################################

#########################################################

#########################################################

#fundamental discriminant finder that Jake came up with, much easier to comprehend.

def fundDisc(N):
    sign = np.sign(N)
    N = abs(N)
    #strip out squares
    for x in range(2, int(math.floor(N ** 0.5)) + 1):
        if 0 == (N % (x**2)):
            return fundDisc(sign * N/(x**2))
    #re-sign the integer
    temp = sign * N
    if 1 == temp % 4:
        return temp
    return 4 * temp

def isFundDisc(N):
    if N == fundDisc(N):
        return True
    else:
        return False


#########################################################

#quick integer sqrt finder, honestly not sure how it works but
#focused on other things right now, will come back once the important
#code is squared away

def isqrt(n):
    x = n
    y = (x + 1) / 2
    while y < x:
        x = y
        y = (x + n / x) / 2
    return int(x)


#########################################################

# kronecker symbol calculator, algorithm taken from Cohen's computational number theory book
# this is already in sage, we should substitute it at some point
def kronecker(a, b):
    if 0 == b:
        if abs(a) != 1:
            return 0
        else:
            return 1
    #remove 2s
    if 0 == a % 2 and 0 == b % 2:
        return 0
    v = 0
    while 0 == b % 2:
        b /= 2
        v += 1
    if 0 == v % 2:
        k = 1
    else:
        if 1==a%8 or 7==a%8:
            k = 1
        elif 3==a%8 or 5==a%8:
            k=-1
    if b < 0:
        b *= -1
        if a < 0:
            k *= -1
    #finished?
    return kronecker2(a,b,k)

def kronecker2(a,b,k):
    if 0 == a:
        if 1 < b:
            return 0
        if 1 == b:
            return k
    v = 0
    while 0 == a % 2:
        v += 1
        a /= 2
    if 1 == v % 2:
        if 1==b%8 or 7==b%8:
            k *= 1
        elif 3==b%8 or 5==b%8:
            k*=-1
    #apply reciprocity
    if 3 == a % 4 and 3 == b % 4:
        k *= (-1)
    r = int(abs(a))
    a = b % r
    b = r
    return kronecker2(a,b,k)


#########################################################


#this function returns divisors of n that take the form p^n, where p is a prime and n is one less than the maximal power of p dividing n

def findPrimes(n):
    #create our set which will contain the prime powers of n with power greater than 1
    #Note: we cannot create a empty set, otherwise it becomes a dictionary, so we start with 1 here and remove it later
    #itereate through integers between 2 and half of n
    primePowers = {1}
    nAb = abs(n)
    for i in range(2, int((nAb/2)+1)):
        #check if i is prime and divides n
        if isPrime(i) and nAb % i == 0:
            #print(i)
            #find largest power (greater than 1) of i that divides n
            y = 1
            while nAb % i ** (y + 1) == 0:
                y += 1
            #when i ** (y+1) does not divide we add i ** y, the greatest power of i dividing n, to our prime powers list
            if y % 2 == 0:
                primePowers.add(i ** y)
            if y % 2 == 1:
                primePowers.add(i ** (y-1))
    #here we remove 1, which was added earlier due to syntax restraints
    primePowers.remove(1)
    return primePowers


#########################################################


# getFundDisc() finds D0 a fundamental discriminant where D=D0f^2
# Jake made a better function. This is obsolete. 

def getFundDisc(D):
    #the way i have constructed the algorithm makes -4 a weird case, so we deal with it seperately for now
    if D == -4:
        return -4
    #here we define D_0, the fundamental discriminant in n
    D_0 = D
    #we iterate through the prime powers given by findPrimes
    for primep in findPrimes(D):
        #we divide out these prime powers to get a squarefree divisor of n
        D_0 = D_0 / (primep)
    #we check the fundamental discriminant condition concerning modulo n (see any standard definition of a fundamental disc)
    #and adjust accordingly
    if D_0 % 4 == 1:
        return D_0
    elif D % 4 == 0:
        return D_0 * 4


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
                MOB = mobius(d)
                KRON = kronecker(D0, d)
                POW = d ** -s
                SIG = sigmaDivisor(1 - (2 * s), int(f / d))
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
def makeCoeff(a,b,c,k):

    if (a,b,c) == (0,0,0):
        return 1

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
            for y in giveReps(-D):
                a = y[0]
                b = y[1]
                c = y[2]
                B[str(y)] = str(makeCoeff(a,b,c,k))
            A[str(-D)] = B
    return A

#########################

def BQFsplitter(a,b,c):
    D = list()
    for n in range(0,a+1):
        a1 = a - n
        for m in range(0,c+1):
            c1 = c - m
            ran = math.floor(2*math.sqrt(a1 * c1))
            for k in range(-int(ran),int(ran)+1):
                S1cand = (a1,k,c1)
                b2 = b-k
                S2cand = (n,b2,m)
                if (b2)**2 <= 4*m*n:
                    D.append((S1cand,S2cand))
    return D


#########################

def bracket(A,S):
    brkt = np.matmul(S.transpose(), np.matmul(A,S))
    return brkt

def matToBQF(A):
    if A[0,1] != A[1,0]:
        pass
        #raise error 
    a = int(A[0,0])
    b = int(2*A[0,1])
    c = int(A[1,1])
    return (a,b,c)

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

        
    def coeff(self,a,b,c):
        if self.type == 'eisen':
            return makeCoeff(a,b,c, self.weight)
        if self.type == 'mult':
            coeff = 0
            F1 = self.constituents[0]
            F2 = self.constituents[1]
            for Q in BQFsplitter(a,b,c):
                Q1 = Q[0]
                Q2 = Q[1]
                coeff += (F1.coeff(Q1[0],Q1[1],Q1[2]) * F2.coeff(Q2[0],Q2[1],Q2[2]))
            return coeff    
        if self.type == 'add':
            F1 = self.constituents[0]
            F2 = self.constituents[1]
            coeff = F1.coeff(a,b,c) + F2.coeff(a,b,c)
            return coeff
        if self.type == 'scale':
            F1 = self.constituents[0]
            c1 = self.constituents[1]
            coeff = c1 * F1.coeff(a,b,c)
            return coeff
        if self.type == 'twist':
            F1 = self.constituents[0]
            p = self.constituents[1]
            #case (i)
            if ((not b % p) == 0 and a % p**4 == 0):
                P_1 = p**(1 - F1.weight) * leg(b,p)
                P_2 = 0
                A=np.matrix([[a,b/2],[b/2,c]])
                for r in mresidues(p, isprime = True):
                    S = np.matrix([[1,Fraction(-b,p)],[0,p]])
                    #print('mark1')
                    newS = matToBQF(bracket(A,S))
                    #print('mark2')
                    P_2 += leg(r,p) * F1.coeff(newS[0],newS[1],newS[2])
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