import numpy as np
import os
from subprocess import PIPE
from subprocess import Popen
from subprocess import run
from typing import Callable  # NOQA

y_best = float('inf')
x_best = np.zeros(0)
x_list = []
y_list = []


def solve(fun, lbounds, ubounds, budget, verbose=True):
    # type: (Callable, np.ndarray, np.ndarray, int, bool) -> None

    global x_best, y_best, x_list, y_list
    y_best = float('inf')
    x_best = np.zeros(0)
    base = os.path.dirname(os.path.abspath(__file__))

    cpp_path = os.path.normpath(os.path.join(base, 'solvers.cpp'))
    exec_path = os.path.normpath(os.path.join(base, 'a.out'))

    run(['g++', cpp_path, '-std=c++11', '-Ofast', '-o', exec_path])

    p = Popen([exec_path], shell=True, stdout=PIPE, stdin=PIPE)
    D = lbounds.shape[0]
    p.stdin.write(bytes(str(D) + ' ' + str(budget) + '\n', 'UTF-8'))
    p.stdin.flush()
    p.stdin.write(bytes(' '.join(map(str, lbounds.tolist())) + '\n', 'UTF-8'))
    p.stdin.flush()
    p.stdin.write(bytes(' '.join(map(str, ubounds.tolist())) + '\n', 'UTF-8'))
    p.stdin.flush()

    for t in range(budget):
        X = np.array(
            p.stdout.readline().decode('utf-8').strip('\n').split(' ')).astype(
                np.float64)
        Y = fun(X)
        p.stdin.write(bytes(str(Y) + '\n', 'UTF-8'))
        p.stdin.flush()
        x_list.append(X)
        y_list.append(Y)
        if y_best > Y:
            y_best = Y
            x_best = X
            if verbose:
                print("iteration:" + str(t) + ",", "X=" + str(X) + ",",
                      "Y=" + str(Y))
