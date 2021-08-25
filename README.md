# Expectation-Maximization-Entropy method

This repository provides the software that goes with the publication [Accurate Detection of Arbitrary Photon Statistic](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.123.153604).

We present a novel algorithm for finding positive solutions to linear inverse problems (LINPOS).
The algorithm is based on expectation-maximization (EM) algorithm [Dempster1977,Vardi1993] weakly [regularized](https://en.wikipedia.org/wiki/Regularization_(mathematics)) by maximum-entropy principle.
The optimal value of the renormalization parameter depends on the amount of data similarly as for other regularization techniques [Caponnetto2006].
Alternatively, the value of the renormalization parameter can be determined adaptively.

We demonstrate the performance of the method in photon statistics retrieval from click statistics data [Hloušek2019].
This problem arises when a state of light is detected by a multi-channel detector with a limited number M of channels.
Each channel is sensed by a binary detector capable of distinguishing two cases - no photons and some photons.
Also, these binary detectors can have a non-unity detection efficiency.
The click statistics (also termed coincidence statistics) c_m, m=0,...,M,
is determined by the input photon statistics c_n, n=1,...\infty,
and the response matrix T_mn of the multi-channel detector, c_m = \sum_n T_mn p_n.
The inverse problem is [ill-posed](https://en.wikipedia.org/wiki/Well-posed_problem) and can be approached by various methods, such as direct pseudoinverse, maximum likelihood estimation, or EME approach. The performance of these methods is compared in Supplemental material of [Hloušek2019].

[Dempster1977]  A. P. Dempster, N. M. Laird, and D. B. Rubin, J. Royal Stat. Soc. B39, 1 (1977).  
[Vardi1993]  Y. Vardi and D. Lee, J. Royal Stat. Soc. B55, 569 (1993).  
[Caponnetto2006] A. Caponnetto and E. De Vito, Found. Comput. Math. 7, 331 (2007).  
[Hloušek2019] J. Hloušek, M. Dudka, I. Straka, and M. Ježek, Phys. Rev. Lett. 123, 153604 (2019). [10.1103/PhysRevLett.123.153604](https://doi.org/10.1103/PhysRevLett.123.153604), preprint: [arXiv:1812.02262](https://arxiv.org/abs/1812.02262).

If you have any questions or feedback, do not hesitate to contact me: Miroslav Ježek, jezek(at)optics.upol.cz, [Quantum Optics Lab Olomouc](http://quantum.opticsolomouc.org/).

## EME algorithm for photon statistics retrieval

Python implementation requirements:
- Numpy
- Scipy (functions binom and factorial)
- Numba (optional)

## eme.py

Minimal example of photon statistics retrieval from measured multichannel coincidences using the EME method.
