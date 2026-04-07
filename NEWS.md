# bigMap 4.7.0, Apr. 2026

## Initial default configuration change

1. Set normalize <- F in src/bdm/glbl.R. By default, bigMap was performing a centering and scaling of the data (see R/bdm_exprt.R and src/zBeta.cpp). This was preventing a correct computation of cosine similarity when passing normalized vector data. Herein, the pre-processing of data is definitely and absolutely left to the user.


# bigMap 4.6.2, Apr. 2024

## Fixes

1. TSNE::getForces() was declared as double when it should be void. This is not detected as an error at compilation time and was not throwing any error until some recent update of R (or gcc++ !!??), when it suddenly starts throwing strange error messages ....

2. Solved bug in bdm.optk.s2nr(). When input data is NULL this function returns an error message string. With ret.optk = TRUE, and using the input bdm instance to store the output, this error message would overwrite the bdm instance, consequently destroying it. Now it throws the error message but still returns the bdm instance.

3. Solved bug in get.levels() function in bdm_auxf.R. bdm.pMap() was failing when Betas are extremely high. Seting the last break resulted in a repeated break value because I was setting it by just adding 1 to the max value of the data. This is solved by setting it to Inf.

4. Improved bdm.qMap()

# bigMap 4.6.1, Feb. 2024

## Fixes

1. Solved bug in s2nr.value() in file R/bdm_s2nr.R.


# bigMap 5.2.2, Apr. 2021

Copy of 4.2.2 implementing the bug found in FitSNE: use distances instead of squared distances in HD.

# bigMap 4.2.2, Mar. 2021

Just a few differences with respect to 4.2.1.

1. Adjusted scheme for perplexity:

	- base_nnS set to bdm.nnSize(base_ppx = 0.03, xppx = 3) (9 neighbors out of 100)
	- base_nX set to 2**20 = 1048576
	- base_epochs set to 2**7 = 128 (number of epochs for nX = base_nX, ppx = base_ppx)
	- epochs <- floor(128 * log(nX**2) /log(base_nX**2) * log(base_nnS) /log(nnS)) = floor(128 * log(nX) /20 * log)
	- base_zSize set to 2**14 = 16384
	- iters <- ceiling(128 * log(zSize**2) /log(base_zSize**2)) = ceiling(128 * log(zSize)/14)

2. Changed learning rate: (dropped factor 1/2 in range(Y))

	- eta = log(zSize -1) * range(Y)

3. Dropped parameter boost. Does not make sense anymore after setting this scheme.

4. Dropped parameter nnExact.

# bigMap 3.9.5, Jan. 2021

The neighbors' set size limit is to strict. Datapoints with high precision form a cloud around the embedding not getting mapped to their right positions. I assume this is because many of their neighbours have Pij = 0 and too few play attractive forces on them and as a result the repulsive forces push them outside of the embedding. This is particularly clear with low perplexities, when the neighbors' sets are smaller. Thus, removed the nSS limit.

Also, split affinity functions into X2P_exact and X2p_apprx. hoping to improve running times.

# bigMap 3.9.4, Dec. 2020

Reimplemented neighbours' set control by means of a distance quantile and a neighbor's set size limit:

	- distance quantile is determined by a xSigma parameter such that:
	Lmax = 1/Bi * log(ppx /Zi) + xSigma**2 /(2 * Bi)

	- xSigma = 0 is used to fix Lmax = L_3ppx

	- the neighbor' set size limit is proportional to 3ppx for the size of the thread, i.e.:
	nSS = max(15, min(3 * ppx * z/n, (z -1)))

The distance quantile is computed in zBeta.cpp() and returned along with Bi in a matrix, not in a list with separated vectors as before. Thus, both parameters are exported together and handled more efficiently.

Added new utility bdm.chknSS() to check the nSS for different values of xSigma.

# bigMap 3.9.3, Dec. 2020

Version 3.9.2 shows NO SIGNIFICANT improvement in running time with respect to 3.9.1 (rather worst). Back to 3.9.1.

Corrected bug in 3.9.1 affmtx.cpp: Li was set to 32/Bi instead of 18/Bi in affMtx::D2P() and affmtx::S2P().

Also, removed theta check in affmtx.cpp so that even for theta <=0 all Lij < 18/Bi yield P[ij] = 0. The counterpart is that W is always computed though not used for theta = 0.

I like the current scheme but testing with P1G it becomes obvious that nearest neighbors' sets up to Lij < 18/Bi results in huge neighbors' sets and unacceptable running times. Added a new parameter xSigma to control the neighbors' set size (defaults to 6, meaning Lij < 6**2/(2*Bi), but allows values down to 3 meaning Lij < 3**2/(2*Bi) = 4.5/Bi).

Using xSigma DOES NOT WORK!!! In P100K, the minimum pair-wise distance is 1. When Bi is higher than xSigma**2/2 we have Lmax = xSigma**2/(2 * Bi) < 1, thus NO Lij is greater than Lmax, and therefore Pj|i = 0 for all j. Again, I'm stepping with the idea that we can not assume anything about the distribution of distances and we can only control the nearest neighbors' sets by fixing a size, or by means of a distance-quantile as I was doing some versions before.

Worth noting that the smallest the neighbors' sets the larger the Pj|i and the stronger the attractive forces, while the intuition is that the precision of the mapping should be lower as we are taking into account less pair-wise interactions. In other words, the effect of using smaller neighbors' sets is equivalent to the effect of using lower values of perplexity relative to the values I was using till now. A side effect is that the algorithm will be faster.

# bigMap 3.9.2, Dec, 2020

Cloned 3.9.1 without recasting of bigMatrix in affmtx.cpp. Just to compare running times with respect to 3.9.1.

# bigMap 3.9.1, Dec, 2020

- Added bdm.restart() with new parameter epochs.

- New chunk&mix scheme (n:dataset size, l:layers, c:threads):

1. epochs = ceiling(sqrt(n * sqrt(l/c)))
2. iters = floor(sqrt(n * sqrt(l/c)))

As previously set (epochs = sqrt(n), iters = sqrt(n * l/c)), iters is decreased with the number of threads but epochs is not. Thus, large n leads to a scheme with too many epochs (too long running time).

Both settings yield the same total amount of iterations:

1. old: sqrt(n) * sqrt(n * l/c) = n * sqrt(l/c)
2. new: sqrt(n * sqrt(l/c)) * sqrt(n * sqrt(l/c)) = n * sqrt(l/c)

The benefit is a reasonable number of epochs for large n. The downside is that given a fix thread size, the number of iterations increases with n to compensate the smaller number of epochs, and might have an effect on the convergence between threads.

Thus, decided to combine both schemes:

1. epochs = ceiling(sqrt(n * sqrt(l/c)))
2. iters = ceiling(sqrt(n * l/c))

This scheme brings both: a reasonable number of epochs for large n, and a constant number of iters for a fixed thread size and, consequently, a smaller number of total iterations, that is:

1. faster running times
2. same degree of convergence among threads as with the previous scheme
3. for small data sets solutions are likely to be not so optimal (unless running more rounds)

# bigMap 3.9.0, Dec, 2020

Computation of Betas with dynamic tolerance and up to Lij < 18/Bi (exp(-18) = 1.522998e-08)

Decided to drop exact computation of affinities and t-SNE. Eventually, keep this option linked to a specific parameter (exact = TRUE, meanwhile I use theta < 0) for very special cases where running time is not an issue and the highest precision might be beneficial (i.e. very small data sets like pBrains).

1. theta < .0   : Pij > 0 for all i,j    ; exact t-SNE
2. theta < .25  : Pij > 0 for Lij < 18/Bi; exact t-SNE
3. theta >= .25 : Pij > 0 for Lij < 18/Bi; aprox t-SNE

# bigMap 3.7.8, Dec, 2020

Changed learning rate to: eta = threads /layers * log(z -1) * 2.

DOES NOT WORK !!!

Changed learning-rate back to eta = range(Y) /2 * log(z -1) * 4, as it was some versions before. Optimizaion is slower but convergence among threads and layers is better.

Betas are computed up to Lij < 16 /Bi and with tolerance log(1 /.99). This leads to significant differences in the accuracy of the betas but not so significant differences in the embedding (see sprsmtx/sierpinki3d/bm378/tol.99 and sprsmtx/sierpinki3d/bm378/tol.99999). But on the other side the computation is much faster and makes it worth.

# bigMap 3.7.7, Dec, 2020

Dynamic computation of betas: betas are computed up to Lij < 32/Bi (as exp(-Bi * 32/Bi) = exp(-32) = 1.266417e-14, almost the .Machine$double.eps = 2.220446e-16).

Abandoned the theta-quantiles: affinities also computed up to Lij < 32/Bi. Nevertheless, still keep the exact computation of affinities (theta == 0) for small datasets.

Changed constructor in affmtx.cpp. Constructor gets SEXPs to X and B and sets them as attributes so that affMtx.methods() can skip this step.

# bigMap 3.7.6, Nov, 2020

Assumed that affinities with Lij < Li(theta) will be very close to zero and can be discarded for the exact version (theta = .0 and Li = 6*ppx). So, removed the computation for theta = .0 in affMtx.cpp. A side effect is that the computation of affinities should be faster as I'm not checking the value of theta for all n*n(-1) pair-wise affinities (could be significant for large thread sizes (??)).

After some thinking about this, my conclusion is that it is a bad idea to discard the possibility of performing the most exact computation possible, because this can be really significant with small datasets (e.g. primate brains). So, left as it was.

Also, after checking with sierpinski3d I realize that it is indeed NOT true that affinities for Lij < Li(6*ppx) are close to zero. In this case, it may be due to the fact that distances (given as short-path-distances) are highly discrete and we might be discarding a large set of affinities for which Lij is exactly Li(6*ppx). So, changed this check to Lij <= Li.

Results with mnist optical digits confirm that adjustment of epochs as of 3.7.4 is not enough for this dataset, so readjusted epochs as:
	epochs <- ceiling(ifelse(N<75000, N**0.50, N**0.45) /bdm$ptsne$boost)
	iters <- ceiling(sqrt(thrd.size) x bdm$ptsne$boost)

Corrected BUG: When starting an 'MPI' cluster check that (length(cl) == threads), as the whole 'MPI' implementation is based on the assumption that we have 1 core (thread 0) running the master process and at least as much cores as to run one working-thread each one (NO MULTI-THREADING).

# bigMap 3.7.4, Nov, 2020

Recasting SEXP to C++ standard pointers as in version 4.0.0, 4.1.0 does not lead to any improvement in running times (amazingly, it rather seems to worsen the performance!), so it does not justify at all the cost of transposing the input data matrix and adapting all the package functionality affected by this transposing. So, I definitely decided to leave this option aside and work with SEXP.

Summary of this version:

	- betas are computed up to 3*ppx
	- default theta is set to 0.5
	- the theta quantile (L_theta)is set to (3*ppx + 2*ppx (1-theta))
	- exact affinities (theta = 0) include all pair-wise distances by thread
	- aprox. affinities (theta >0.2) include pair-wise distances such that Lij < L(thtq)
	- affinities are normalized by thread
	- learning rate is set to 2*(log(z) + log(z-1))
	- allows computing the output for multiple perplexities without stopping the cluster (and reloading data each time) and computing the k-ary neighbourhood preservation (setting qlty = TRUE) for comparison purposes
	- epochs/iters are adjusted for large datasets as:
	epochs <- ceiling(ifelse(N<10000, N**0.50, N**0.45) /bdm$ptsne$boost)
	iters <- ceiling(sqrt(thrd.size) x bdm$ptsne$boost)
	- epochs can be forced (even for multiple perplexities)

Up to this points, results confirm that output noise is only due to the effect of running multiple partial tSNEs. Increasing the thread-size (either using less threads or more layers) reduces the noise in the output and with thread-size == dataset-size the output is similar to openTSNE implementations (even better I would say).

# bigMap 3.7.2. Oct, 2020

Going back to thread normalization of P and Q.

# bigMap 3.7.1. Oct, 2020

Previous version resulted not so fast and not so accurate as expected ... Instead of using a fix set of nearest neighbours, given by thtQ = 3pxx independently of theta (following V.d.Maaten), I use
thtQ = ppx + 2ppx (1 -theta).

Number of epochs changed from sqrt(N) to N**0.45.

# bigMap 3.7.0. Set, 2020

Changes in the computation of affinities. When computing betas, row normalization factors are saved. Afterwards, thread affinities are straightforwardly normalized using these factors. Thus, computation of affinities is faster and do not depend on the particular set of observations on the thread-set. In other words, distributions P and Q are GLOBAL at each thread, though we consider only chunks of them. But, how can we manage a GLOBAL version of Q at thread level ??

Instead of bootstrapping the data set and making independent runs of the tSNE, we make partial tSNEs where affinities are always relative to the global P and Q distributions.

This change induces a major change in the structure of code, where computation of affinities is detached from tSNE. Also, following Van der Maaten, input affinities are computed only for the set of neighbours in the quantile 3pxx, even in the case of the exact algorithm. The result is a faster and more simple code.

# bigMap 3.6.0. Set, 2020

Added support for sparse matrices.

# bigMap 3.5.0. July, 2020

Raw data matrix is dropped from the bdm object. Raw data matrix is kept aside as a bigmemory::big.matrix() object by attaching its reference to the bdm instance. This saves memory space and a significant amount of computation time when passing a bdm instance from one function to the other.

# bigMap 3.3.2. June, 2020

Quality measures based on Hellinger divergence.

# bigMap 3.3.0. May, 2020

Index of structure preservation, new implementation.

# bigMap 3.2.0. Gen, 2020

Implementing ptSNE tree-based approximations: beta-quantiles (hd-space), Barnes-Hut quad-tree (ld-space)

# bigMap 3.1.0. Dec, 2019

Reimplementation of the computation of the optimal number of input dimensions.
Changes in PCA data preprocessing (bdm_data.R)

# bigMap 3.0.0. Dec, 2019.

Faster computation of the similarity matrix.
New feature: Index of structure preservation.

# bigMap 2.3.0. Nov, 2019.

New feature: computation of optimal number of input-dimensions

# bigMap 2.2.0. Oct, 2019.

New feature: S2NR-based merging of clusters.
Improved quantile-Map and density-Map plots.
The generic bdm.plot() has been replaced by particular plot functions bdm.ptsne.plot(), bdm.pakde.plot(9 and bdm.wtt.plot().

# bigMap 2.1.0. Feb, 2019.

Solved for staged install.

# bigMap 2.0.0. Jan, 2019.

First submission.
