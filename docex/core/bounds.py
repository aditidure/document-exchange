"""
bounds.py -- the theoretical quantities the measurements are compared with (report, Sec. II).

    theorem1(n, t)            4t log n                       leading term of Theorem 1
    theorem2(n, t)            (4t - 1) log n                 leading term of Theorem 2
    lemma3_log_ball(n, t, k)  log2 of Eq. (2): prod_{i=0}^{2t}(n + ik) k^{2t} 3^{2kt}  (upper bound on |C_t(x)|)
    estimated_tag(n, t, k)    2 log2 of that bound: the tag length Theorem 1 guarantees, used where a string is too
                              long to encode by enumeration (whole real files in the extension)
"""
import math


def theorem1(n, t):
    return 4 * t * math.log2(n)


def theorem2(n, t):
    return (4 * t - 1) * math.log2(n)


def lemma3_log_ball(n, t, k):
    return sum(math.log2(n + i * k) for i in range(2 * t + 1)) + 2 * t * math.log2(k) + 2 * k * t * math.log2(3)


def estimated_tag(n, t, k):
    return 2 * lemma3_log_ball(n, t, k)
