"""Compute the Interpolation Polynomial"""
import numpy as np

def PolValue(coeffs, t):
    '''
    Value of interpolation polynomial in point t
        degree       degree of polynomial
        nodesX  Nodes of interpolation: x0, x1, ..., xn
        coeffs  Coefficients of polynomial
        t       Point of interpolation
    '''
    s = 0
    for i in range(0, len(coeffs)):
        s += (coeffs[-(i+1)]) * t**i;
    return s

def computePol(nodes,degree):
    '''
    Compute the coefficients of
    interpolation polynomial
        n       Number of nodes - 1 == degree of polynomial
        nodes   Nodes of interpolation:
                [(x0, y0), ..., (xn, yn)]
        Return:
        coeffs  Coefficients of Newton polynomial
    '''
    Y=np.array([p[1] for p in nodes])
    coeffs = [0.0]*(degree + 1)
    A = np.ones((len(coeffs),len(nodes)))
    for i in range(A.shape[0]):
        A[i]=[p[0]**(degree-i) for p in nodes]
    A = A.T
    print(A)
    pinvA=np.dot(np.linalg.inv(np.dot(A.T,A)), A.T)
    print(pinvA)
    coeffs_T=np.dot(pinvA, Y.T)
    return coeffs_T