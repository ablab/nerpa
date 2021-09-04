import numpy as np

from typing import List


def my_log(a):
    if a == 0:
        return -np.inf
    return np.log(a)


def my_exp(a):
    if a == -np.inf:
        return 0
    res = np.exp(a)
    if isinstance(res, np.ndarray):
        return res[0]
    return res


def log_add_exp(l: List[float]):
    l = list(filter(lambda x: x != -np.inf, l))
    if len(l) == 0:
        return -np.inf
    if len(l) == 1:
        return l[0]
    l = sorted(l)
    res = np.logaddexp(l[0], l[1])
    for i in l[2:]:
        res = np.logaddexp(res, i)
    return res