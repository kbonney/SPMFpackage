load('spmf.sage')

X=BinaryQF(81,1,1)
Y=BinaryQF(1,0,1)
Z=BinaryQF(0,0,0)

A = matrix([[81,1/2],[1/2,1]])
S = matrix([[1,3],[0,-1/3]])

T1=multiplyForms(SPMF(10),SPMF(10))
T2=multiplyForms(SPMF(8),SPMF(12))
T3=SPMF(20)
C=addForms(T1,T2)
D=twistForms(C,3)


