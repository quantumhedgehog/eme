# Expectation-Maximization-Entropy method

We present a novel algorithm for finding positive solutions to linear inverse problems (LINPOS).
The algorithm is based on expectation-maximization (EM) algorithm [] weakly regularized by maximum-entropy principle.
The optimal value of the renormalization parameter depends on the amount of data similarly as for other regularization techniques [Caponnetto2006].
Alternatively, the value of the renosrmalization parameter can be determined adaptively.

We demonstrate the performance of the method in photon statistics retrieval from click statistics data [Hloušek2019].
This inverse problem arises when a state of light is detected by a multi-channel detector with a limited number M of channels.
Each channel is sensed by a binary detector capable of distinguishing two cases - no photons and some photons.
Also, these binary detectors can have a non-unity detection efficiency.
The click statistics (also termed coincidence statistics) c_m, m=0,...,M,
is determined by the input photon statistics c_n, n=1,...\infty,
and the response matrix T_mn of the multi-channel detector, c_m = \sum_mn T_mn p_n.
The inverse problem is [ill-posed](https://en.wikipedia.org/wiki/Well-posed_problem) and can be approached by several methods,
see supplemental material [Hloušek2019].
