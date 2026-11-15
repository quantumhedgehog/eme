# cython: boundscheck=False, wraparound=False, cdivision=True
import numpy as np
from libc.math cimport logl, sqrtl, isnan

def EME(long double[:, ::1] det, long double[::1] c,
        long double[::1] lval, long double[::1] epsval,
        long long iterations=100000000):
    cdef Py_ssize_t mMax = det.shape[0], nMax = det.shape[1]
    cdef Py_ssize_t m, n
    cdef long long iteration
    cdef long double l = lval[0], epsilon = epsval[0]
    cdef long double predicted, entropy, em, ent, norm, dist, delta
    cdef bint entropy_nan
    cdef long double[::1] pn = np.empty(nMax, dtype=np.longdouble)
    cdef long double[::1] nxt = np.empty(nMax, dtype=np.longdouble)
    cdef long double[::1] aux = np.empty(mMax, dtype=np.longdouble)
    cdef long double[::1] back = np.empty(nMax, dtype=np.longdouble)

    for n in range(nMax):
        pn[n] = (<long double>1)/nMax

    for iteration in range(iterations):
        for m in range(mMax):
            predicted = 0
            for n in range(nMax):
                predicted += det[m,n]*pn[n]
            if predicted == 0 and c[m] == 0:
                aux[m] = 0
            else:
                aux[m] = c[m] / predicted

        for n in range(nMax):
            back[n] = 0
        for m in range(mMax):
            for n in range(nMax):
                back[n] += det[m,n]*aux[m]

        entropy = 0
        entropy_nan = False
        if l != 0:
            for n in range(nMax):
                if pn[n] <= 0:
                    entropy_nan = True
                    break
                entropy += pn[n]*logl(pn[n])

        norm = 0
        for n in range(nMax):
            em = back[n] * pn[n]
            if l == 0 or entropy_nan:
                ent = 0
            else:
                ent = l * (logl(pn[n]) - entropy) * pn[n]
            nxt[n] = em - ent
            norm += nxt[n]

        dist = 0
        for n in range(nMax):
            nxt[n] /= norm
            delta = nxt[n] - pn[n]
            dist += delta*delta
        if sqrtl(dist) <= epsilon:
            return np.asarray(nxt).copy(), iteration + 1

        for n in range(nMax):
            pn[n] = nxt[n]

    return np.asarray(pn).copy(), iterations
