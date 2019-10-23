import numpy as np
import solver.Limited_GP as gp


def objective(x):
    x_optimal = np.array([-1, 3])
    y_optimal = 37
    return y_optimal + np.linalg.norm(x - x_optimal, ord=2)


lbound = np.full(2, -5)
ubound = np.full(2,  5)
budget = 100

gp.solve(objective, lbound, ubound, budget)
print("min:\t", gp.y_best)
print("argmin:\t", gp.x_best)
