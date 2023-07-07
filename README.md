# SPMFpackage
This repository houses a collection of code designed to implement the twisting algorithm proposed in [Twisting of Siegel paramodular forms](https://arxiv.org/abs/1404.4596).
The core class `SPMF` represents Siegel paramodular forms. Standard vector space operations are implemented for SPMF as well as a method for twisting.
This repository is only in an alpha stage, so many things are unoptimized or poorly documented. I no longer work on this project, but I can help answer any questions about what is here.


# fundDisc(N) and isFundDisc(N)
fundDisc(N) takes an integer N and returns a fundamental discriminant dividing it. This is important for calculating Cohen's class number function which is involved in the formula for coefficients of the Siegel Eisenstein series. 

isFundDisc(N) simply tells us if N is already a fundamental discriminant.

## kronecker(a,b) 
kronekcer(a,b) calculates the kronecker symbol with a over b. This is used in the formula for coefficients of the Siegel Eisenstein series.

## makeCoeff(a,b,c,k)
makeCoeff(a,b,c,k) takes a binary quadratic form (a,b,c) and a weight k and returns the coefficient  of the Siegel Eisenstein series of weight k indexed by (a,b,c). The following code block gives an example of the calculation of the coefficient indexed by (1,1,1) for the Siegel Eisenstein series of weight 4.
```
┌────────────────────────────────────────────────────────────────────┐
│ SageMath version 8.8, Release Date: 2019-06-26                     |
│ Using Python 2.7.15. Type "help()" for help.                       |
└────────────────────────────────────────────────────────────────────┘
sage: import functionPile as fp
sage: fp.makeCoeff(1,1,1,4)
13440
sage: 
```

### genJson(k,N) 
genJson(k,N) takes a weight k and a positive integer N, and returns an indexed list of coefficients for the Siegel Eisenstein series of weight k. The integer N is a cap on the discrimiant of the binary quadratic forms which will be listed. Currently only the coefficients for primitive forms are printed. In practice, this function is not used directly out of the module, but instead is best used via *makeEisenJson.sage*. A usage might look like:
```
kirk@silver:~/twists/sage$ sage  makeEisenJson.sage 4 12
{
 "0": {}, 
 "-3": {
  "x^2 + x*y + y^2": "13440.0"
 }, 
 "-4": {
  "x^2 + y^2": "30240.0"
 }, 
 "-7": {
  "x^2 + x*y + 2*y^2": "138240.0"
 }, 
 "-8": {
  "x^2 + 2*y^2": "181440.0"
 }, 
 "-11": {
  "x^2 + x*y + 3*y^2": "362880.0"
 }
}

```

This gives us the coefficients for the weight k form indexed by primitive forms of fundamental discriminant up through -12. 

### BQFsplitter(a,b,c)
BQFsplitter(a,b,c) takes a positive semidefinite binary quadratic form and splits it into all possible partitions S1 + S2 = S where S1 and S2 are also positive semidefinite binary quadratic forms. This is important for the multiplication of Siegel paramodular forms. Here is an example of its usage on the binary quadratic form (2,0,2):
```
┌────────────────────────────────────────────────────────────────────┐
│ SageMath version 8.8, Release Date: 2019-06-26                     |
│ Using Python 2.7.15. Type "help()" for help.                       |
└────────────────────────────────────────────────────────────────────┘
sage: import functionPile as fp
sage: fp.BQFsplitter(2,0,2)
[((2, 0, 2), (0, 0, 0)),
 ((2, 0, 1), (0, 0, 1)),
 ((2, 0, 0), (0, 0, 2)),
 ((1, 0, 2), (1, 0, 0)),
 ((1, -2, 1), (1, 2, 1)),
 ((1, -1, 1), (1, 1, 1)),
 ((1, 0, 1), (1, 0, 1)),
 ((1, 1, 1), (1, -1, 1)),
 ((1, 2, 1), (1, -2, 1)),
 ((1, 0, 0), (1, 0, 2)),
 ((0, 0, 2), (2, 0, 0)),
 ((0, 0, 1), (2, 0, 1)),
 ((0, 0, 0), (2, 0, 2))]
sage: 
```

### class SPMF(k)




