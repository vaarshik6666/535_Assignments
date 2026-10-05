import numpy as np
import cvxpy as cp

def calc_Jacobian(x, u, param):

    L_f = param["L_f"]
    L_r = param["L_r"]
    dt   = param["h"]

    psi = x[2]
    v   = x[3]
    delta = u[1]
    a   = u[0]

    # Jacobian of the system dynamics
    A = np.zeros((4, 4))
    B = np.zeros((4, 2))

    #############################################################################
    #                    TODO: Implement your code here                         #
    #############################################################################
    
    constant = (L_r / (L_r + L_f))
    beta = np.arctan(constant * np.arctan(delta))
    diff_beta_delta = constant / ((1 + (delta**2)) * (1 + ((constant ** 2) * (np.arctan(delta) ** 2))))
    
    A[0, 0] = 1
    A[0, 2] = -1 * dt * v * np.sin(psi + beta)
    A[0, 3] = dt * np.cos(psi + beta)
    
    A[1, 1] = 1
    A[1, 2]  = dt * v * np.cos(psi+beta)
    A[1, 3] = dt * np.sin(psi + beta)
    
    A[2, 2] = 1
    A[2, 3] = (dt * np.arctan(delta)) / ((((((L_r**2) * (np.arctan(delta)**2)) / ((L_f + L_r)**2)) + 1) ** 0.5) * (L_f + L_r))
    
    A[3, 3] = 1
    
    B[0, 1] = -1 * (dt * L_r * v * np.sin(psi + np.arctan(constant * np.arctan(delta)))) / (((delta**2) + 1) * (((constant**2) * (np.arctan(delta)**2)) + 1) * (L_f + L_r))
    
    B[1, 1] = (dt * L_r * v * np.cos(psi + np.arctan(constant * np.arctan(delta)))) / (((delta**2) + 1) * (((constant**2) * (np.arctan(delta)**2)) + 1) * (L_f + L_r))
    
    B[2, 1] = (dt * v) / (((delta**2) + 1) * ((((constant**2) * (np.arctan(delta)**2)) + 1)**(3/2)) * (L_f + L_r))
    
    B[3, 0] = dt
    
    #############################################################################
    #                            END OF YOUR CODE                               #
    #############################################################################

    return [A, B]

def LQR_Controller(x_bar, u_bar, x0, param):
    len_state = x_bar.shape[0] # gives no.of rows in x_bar
    len_ctrl  = u_bar.shape[0] # gives no.of rows in u_bar
    dim_state = x_bar.shape[1] # gives no.of columns in x_bar
    dim_ctrl  = u_bar.shape[1] # gives no.of columns in u_bar

    n_u = len_ctrl * dim_ctrl
    n_x = len_state * dim_state
    n_var = n_u + n_x

    n_eq  = dim_state * len_ctrl # dynamics
    n_ieq = dim_ctrl * len_ctrl  # input constraints

    
    #############################################################################
    #                    TODO: Implement your code here                         #
    #############################################################################

    # define the parameters
    Q = np.eye(4)  * 10
    R = np.eye(2)  * 5
    Pt = np.eye(4) * 10

    # define the cost function
    P = np.zeros((124, 124))
    P[0:80, 0:80] = np.kron(np.eye(20), Q)
    P[80:84, 80:84] = Pt
    P[84:124, 84:124] = np.kron(np.eye(20), R)
    q = np.zeros((124, 1))
    
    # define the constraints
    A = np.zeros((n_eq, n_var))
    b = np.zeros(n_eq)
    
    for k in range(len_ctrl):
        Ak, Bk = calc_Jacobian(x_bar[k, :], u_bar[k, :], param)

        A[k*dim_state:(k+1)*dim_state,
        k*dim_state:(k+1)*dim_state] = Ak

        A[k*dim_state:(k+1)*dim_state,
        (k+1)*dim_state:(k+2)*dim_state] = -np.eye(dim_state)

        A[k*dim_state:(k+1)*dim_state,
        n_x + k*dim_ctrl:n_x + (k+1)*dim_ctrl] = Bk
    
    # Define and solve the CVXPY problem.
    x = cp.Variable(n_var)
    objective = cp.Minimize(0.5 * cp.quad_form(x, P) + q.flatten() @ x)
    constraints = [
    A @ x == b,
    x[:dim_state] == x0 - x_bar[0, :]
]
    prob = cp.Problem(objective, constraints)
    prob.solve(verbose=False, max_iter=10000)


    #############################################################################
    #                            END OF YOUR CODE                               #
    #############################################################################

    u_act = x.value[n_x:n_x + dim_ctrl] + u_bar[0, :]
    return u_act

def CMPC_Controller(x_bar, u_bar, x0, param):
    len_state = x_bar.shape[0]
    len_ctrl  = u_bar.shape[0]
    dim_state = x_bar.shape[1]
    dim_ctrl  = u_bar.shape[1]
    
    n_u = len_ctrl * dim_ctrl
    n_x = len_state * dim_state
    n_var = n_u + n_x

    n_eq  = dim_state * len_ctrl # dynamics
    n_ieq = dim_ctrl * len_ctrl # input constraints

    a_limit = param["a_lim"]
    delta_limit = param["delta_lim"]
    
    #############################################################################
    #                    TODO: Implement your code here                         #
    #############################################################################
    
    # define the parameters
    Q = np.eye(4)  * 10
    R = np.eye(2)  * 5
    Pt = np.eye(4) * 10
    
    # define the cost function
    P = np.zeros((n_var, n_var))
    P[0:80, 0:80] = np.kron(np.eye(20), Q)
    P[80:84, 80:84] = Pt
    P[84:124, 84:124] = np.kron(np.eye(20), R)
    q = np.zeros((n_var, 1))
    
    # define the constraints
    A = np.zeros((n_eq, n_var))
    b = np.zeros(n_eq)
    for k in range(len_ctrl):
        Ak, Bk = calc_Jacobian(x_bar[k, :], u_bar[k, :], param)

        A[k*dim_state:(k+1)*dim_state,
        k*dim_state:(k+1)*dim_state] = Ak

        A[k*dim_state:(k+1)*dim_state,
        (k+1)*dim_state:(k+2)*dim_state] = -np.eye(dim_state)

        A[k*dim_state:(k+1)*dim_state,
        n_x + k*dim_ctrl:n_x + (k+1)*dim_ctrl] = Bk
    
    G = np.zeros((n_ieq, n_var))
    ub = np.zeros(n_ieq)
    lb = np.zeros(n_ieq)

    for k in range(len_ctrl):
        G[2*k, n_x + k*dim_ctrl] = 1
        G[2*k+1, n_x + k*dim_ctrl + 1] = 1

        lb[2*k] = -a_limit - u_bar[k, 0]
        ub[2*k] =  a_limit - u_bar[k, 0]

        lb[2*k+1] = -delta_limit - u_bar[k, 1]
        ub[2*k+1] =  delta_limit - u_bar[k, 1]

    # Define and solve the CVXPY problem.
    x = cp.Variable(n_var)
    objective = cp.Minimize( 0.5 * cp.quad_form(x, P) + q.flatten() @ x)
    constraints = [
        A @ x == b,
        x[:dim_state] == x0 - x_bar[0, :],
        G @ x <= ub,
        G @ x >= lb
    ]
    prob = cp.Problem(objective, constraints)
    prob.solve(verbose=False, max_iter=10000)

    #############################################################################
    #                            END OF YOUR CODE                               #
    #############################################################################
    
    u_act = x.value[n_x:n_x + dim_ctrl] + u_bar[0, :]
    return u_act