#  eme_longdouble.py
#  
#  Photon statistics retrieval from measured multichannel coincidences using the EME method.
#  C_mn matrix is computed using mpmath library with 400 digits
#  and EME iterations are evaluated in longdouble with approximately 19 digits of precision

import numpy as np
import mpmath as mp

# response function (matrix) of multichannel detector
# eta - quantum efficiency
# mMax - nuber of detection channels
# nMax - nuber of elements of photon statistics
def DET(mMax, nMax, eta):
    det = np.zeros((mMax + 1, nMax + 1), dtype=np.longdouble)
    with mp.workdps(400):
        eta = mp.mpf(str(eta))
        for m in range(mMax + 1):
            for n in range(nMax + 1):
                if m > n:
                    det[m, n] = 0
                elif m < n:
                    summary = mp.fsum(
                        (-1)**j * mp.binomial(m, j)
                        * (1 - eta + (m - j)*eta/mMax)**n
                        for j in range(m + 1)
                    )
                    det[m, n] = np.longdouble(mp.nstr(mp.binomial(mMax, m)*summary, 40))
                else:
                    x = (eta/mMax)**n * mp.factorial(mMax)/mp.factorial(mMax - n)
                    det[m, n] = np.longdouble(mp.nstr(x, 40))
    return det

# alternatively, for intermediate sizes of C matrix
# and when speed is crucial, recurrent version of the evaluation 
# can be used
def DET_recurrent(mMax, nMax, eta):
    det = np.zeros((mMax + 1, nMax + 1), dtype=np.longdouble)
    det[0, 0] = 1

    m = np.arange(mMax + 1, dtype=np.longdouble)
    q = eta / np.longdouble(mMax)

    stay = 1 - eta + m*q
    grow = (mMax - m[1:] + 1)*q

    for n in range(nMax):
        det[:, n+1] = stay*det[:, n]
        det[1:, n+1] += grow*det[:-1, n]

    return det

# photon statistics retrieval using EME algorithm in longdouble precision
# iterations - the maximum number of iterations (actual number is much lower, based on "epsilon")
# epsilon - the distance between n-th and (n+1)-th iterations when the process is stopped
# c - click statistics
# pn - photon statistics
def EME(mMax, nMax, eta, det, l, c):
    iterations = 10**12
    epsilon = np.longdouble('1e-18')
    pn = np.array([np.longdouble(1)/(nMax + 1)]*(nMax + 1), dtype=np.longdouble)
    iteration = 0
    with np.errstate(divide='ignore', invalid='ignore'):
        while iteration < iterations:
            aux = c / np.dot(det, pn)
            aux[np.isnan(aux)] = 0
            EM = np.dot(aux, det) * pn
            E = l * (np.log(pn) - np.sum(pn * np.log(pn))) * pn
            E[np.isnan(E)] = 0
            EME = EM - E
            EME /= np.sum(EME)
            dist = np.sqrt(np.sum((EME - pn)**2))
            if dist <= epsilon:
                print("number of iterations:", iteration)
                break
            else:
                pn = EME
                iteration += 1
    return EME

# numerical values
mMax = 10
nMax = 50
eta = np.longdouble('0.5')
l = np.longdouble('1e-4')
# regularization paramater lambda should be smaller than
# 1/sqrt(total number of click of PNRD data, i.e. total number state copies measured)
# and (probably) larger than 1/copies
copies = int(1e7)
alpha = np.longdouble('0.75') # alpha=0.5 gives 1/sqrt(copies), alpha=1 gives 1/copies
l = np.longdouble(copies)**(-alpha)
print("lambda:", l)

# evaluation of C_mn matrix
det = DET(mMax, nMax, eta)
#det = DET_recurrent(mMax, nMax, eta)

# measured click statistics
c = np.array([
    6.73794700e-03, 4.37104954e-02, 1.27601677e-01,
    2.20741125e-01, 2.50599060e-01, 1.95082729e-01,
    1.05461930e-01, 3.90945126e-02, 9.51054071e-03,
    1.37104223e-03, 8.89424261e-05
], dtype=np.longdouble)
# alternatively inputted as strings to be on the safe side
# to avoid conversion to double when experimental data are
# also evaluated in longdouble precision
#c = np.array([
#    '6.73794700e-03', '4.37104954e-02', '1.27601677e-01',
#    '2.20741125e-01', '2.50599060e-01', '1.95082729e-01',
#    '1.05461930e-01', '3.90945126e-02', '9.51054071e-03',
#    '1.37104223e-03', '8.89424261e-05'
#], dtype=np.longdouble)
c /= np.sum(c)
# ideally the click statistics should be evaluated including
# normalization directly from PNRD counts in longdouble precision

# retrieval
p = EME(mMax, nMax, eta, det, l, c)

np.set_printoptions(floatmode='unique')
print("reconstructed statistics:")
print(p)

# forward prediction of the click statistics
c_rec = np.dot(det, p)
dc = c_rec - c

print("total variation distance of the input c and PNRD(EME(c)):", np.sum(np.abs(dc))/2)
