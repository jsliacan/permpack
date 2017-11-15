import numpy as np
import timeit
import sys


# 7 part construction weights
a = [0.2174127723536347308692444843,
     0.0170598057899242722740620549,
     0.0516101402487892270230230972,
     0.4340722809873864994312953007,
     0.1479895625950390496250611829,
     0.0764457255805656971383351365,
     0.0554097124446605236389787433]

def phi7(a):
    """
    Compute density of a construction (Walter's) with 7 parts, each
    having weights according to vector a (parts are labelled
    left-to-right).

    """

    # density 1342 with 3,4 in the same part
    N1 = a[0]*a[6]*(a[1]^2+a[2]^2+a[3]^2+a[4]^2+a[5]^2)/2 + a[0]*a[5]*(a[2]^2+a[3]^2+a[4]^2)/2 + a[0]*a[4]*a[3]^2/2 + a[1]*a[5]*(a[2]^2+a[3]^2+a[4]^2)/2 + a[1]*a[4]*a[3]^2/2 + a[2]*a[3]^2/2*a[4]
    # density of 1342 with 3,4 in different parts
    N2 = a[0]*a[6]*a[1]*(a[2]+a[3]+a[4]+a[5]) + a[0]*a[6]*a[2]*(a[3]+a[4]) + a[0]*a[5]*a[2]*(a[3]+a[4]) + a[1]*a[5]*a[2]*(a[3]+a[4])
    # density of 1342 with 342 in the same part
    N3 = a[0]*a[6]^3/6*n(2*sqrt(3)-3,100)

    return 24*(N1+N2+N3)/(1-a[0]^4)


print "p(1342) >", phi7(a)


