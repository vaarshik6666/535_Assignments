import numpy as np

def ACC_Controller(t, x, param):
    vd = param["vd"]
    v0 = param["v0"]
    m = param["m"]
    Cag = param["Cag"]
    Cdg = param["Cdg"] 

    # cost function and constraints for the QP
    P = np.zeros((2,2))
    q = np.zeros([2, 1])
    A = np.zeros([5, 2])
    b = np.zeros([5])
    
    #############################################################################
    #                    TODO: Implement your code here                         #
    #############################################################################

    # set the parameters
    lam = 0.5
    alpha = 0.5
    w = 1000.0

    # construct the cost function
    P[0, 0] = 2
    P[1, 1] = 2 * w
    
# construct the constraints

    D, v = x[0], x[1]

    h = (v - vd)**2 / 2

    dv = max(v - v0, 0.0)

    B = D - (v - v0)**2 / (2*Cdg) - 1.8*v

    A[0, 0] = (v - vd) / m
    A[0, 1] = -1
    A[1, 0] = -(1.8 + (v - v0) / Cdg) / m
    A[2, 0] = 1 / m
    A[3, 0] = -1 / m
    A[4, 1] = -1

    b[0] = -lam * h
    b[1] = (v0 - v) + alpha * B
    b[2] = Cag
    b[3] = Cdg
    b[4] = 0

    #############################################################################
    #                            END OF YOUR CODE                               #
    #############################################################################
    
    return A, b, P, q