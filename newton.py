import math as math
import numpy as np
import numpy.linalg as nla
import scipy.misc as sm
from gauss_pivot import gaussPivot
from input import damping, h_val

def swapRows(v,i,j):
    if len(v.shape) == 1:
        v[i],v[j] = v[j],v[i]
    else:
        v[[i,j],:] = v[[j,i],:]
def swapCols(v,i,j):
    v[:,[i,j]] = v[:,[j,i]]

def jacobian(f,x, h = h_val):
    n = len(x)
    jac = np.zeros((n,n))
    f0 = f(x)
    for i in range(n):
        temp = x[i]
        x[i] = temp + h
        f1 = f(x)
        x[i] = temp
        jac[:,i] = (f1 - f0)/h
    return jac,f0

def newtonRaphson2(f,x,tol=1.0e-3, full_output = True, damping_val = damping):
    for i in range(500):
        jac,f0 = jacobian(f,x)
        if np.sqrt(np.dot(f0,f0)/len(x)) < tol: return x
        dx = gaussPivot(jac,-f0)
        #dx = nla.solve(jac, -f0)
        x = x + damping_val * dx
        if (full_output):
            print("iteration ", i, "value = ", x)
        if np.sqrt(np.dot(dx,dx)) < tol*max(max(abs(x)),1.0):
            return x
    print("Too many iterations")
