import numpy as np
from scipy.optimize import minimize

# Define the objective function
def objective(x):
    a1, a2, a3, ap1, ap2, ap3, b1, b2, b3, bp1, bp2, bp3 = x
    return -a1*b1 - a2*b2 - a3*b3 - ap1*b1 - ap2*b2 - ap3*b3 + a1*bp1 + a2*bp2 + a3*bp3 - ap1*bp1 - ap2*bp2 - ap3*bp3

# Define constraints
def constraint1(x):
    return np.sum(x[0:3]**2) - 1  # a_1^2 + a_2^2 + a_3^2 = 1

def constraint2(x):
    return np.sum(x[3:6]**2) - 1  # a'_1^2 + a'_2^2 + a'_3^2 = 1

def constraint3(x):
    return np.sum(x[6:9]**2) - 1  # b_1^2 + b_2^2 + b_3^2 = 1

def constraint4(x):
    return np.sum(x[9:12]**2) - 1  # b'_1^2 + b'_2^2 + b'_3^2 = 1

# Constraints imposing orthogonality between vectors on the same detector
# def constraint_orth_as(x):
#         return np.dot(x[0:3],x[3:6])
    
# def constraint_orth_bs(x):
#     return np.dot(x[9:12],x[6:9])

# Initial guess (random unit vectors)
x0 = np.random.rand(12)
x0[:3] /= np.linalg.norm(x0[:3])
x0[3:6] /= np.linalg.norm(x0[3:6])
x0[6:9] /= np.linalg.norm(x0[6:9])
x0[9:12] /= np.linalg.norm(x0[9:12])

# Constraints dictionary
constraints = [
    {'type': 'eq', 'fun': constraint1},
    {'type': 'eq', 'fun': constraint2},
    {'type': 'eq', 'fun': constraint3},
    {'type': 'eq', 'fun': constraint4},
    # {"type": "eq", "fun": constraint_orth_bs},
    # {"type": "eq", "fun": constraint_orth_as}
]

# Solve the optimization problem
result = minimize(objective, x0, constraints=constraints, method='SLSQP')

# Output result
print("Optimal solution:", result.x)
print("Optimal value:", -result.fun)

# Calculate the angles between each pair of vectors from "a dot b = |a||b|cos(theta)""
theta_ab = np.arccos(result.x[0]*result.x[6] + result.x[1]*result.x[7] + result.x[2]*result.x[8])
theta_abp = np.arccos(result.x[0]*result.x[9] + result.x[1]*result.x[10] + result.x[2]*result.x[11])
theta_apb = np.arccos(result.x[3]*result.x[6] + result.x[4]*result.x[7] + result.x[5]*result.x[8])
theta_apbp = np.arccos(result.x[3]*result.x[9] + result.x[4]*result.x[10] + result.x[5]*result.x[11])
theta_aap = np.arccos(result.x[0]*result.x[3] + result.x[1]*result.x[4] + result.x[2]*result.x[5])
theta_bbp = np.arccos(result.x[6]*result.x[9] + result.x[7]*result.x[10] + result.x[8]*result.x[11])

# Prints angles in fractions of pi
angles = [theta_ab/np.pi, theta_apb/np.pi, theta_abp/np.pi, theta_apbp/np.pi, theta_aap/np.pi, theta_bbp/np.pi]
print(angles)
