L-PAC: MATLAB function for determining absolute copy number alterations (CNAs) from multi-region data
L-PAC leverages information from high-purity samples to improve the
inference of CNAs from low-purity samples from the same tumor.
Please note that the description below is for default L-PAC corresponding to the MATLAB function  "LPAC.m". 
A generalized version of which gives the user more flexibility is provided as the MATLAB function "LPAC_generalized.m". The description of how to use the generalized version of LPAC can be found in the MATLAB file "LPAC_generalized.m".

PROCEDURE
1) Single-sample inference of absolute CNAs by the CNH.m algorithm (available at https://github.com/dmmiedema/CNH)
2) Filtering of samples for which single-sample inference in step (1) failed based on the output purity == 1 and ploidy < 2.5.
3) Identification of the best sample from the multi-region data of a tumor as the sample with the lowest intratumor heterogeneity as defined by CNH in step (1), and for which inference was succesfull according to the criteria in step (2).
4) Alignment of relative CNAs to the absolute CNAs in a minimization procedure analogous to that described in Van Dijk et al., Nature Communications 2021, for CNH. But with minimization of candidate absolute CNA values to the best sample identified in step (3), instead of to integer values as in step (1).

DEPENDENCY

1) CNH.m           available at https://github.com/dmmiedema/CNH. Description in Van Dijk et al., Nature Communications 2021 (https://doi.org/10.1038/s41467-021-23384-6).

INPUT

1) bin_val      Nb x Ns size array containing relative CNA values of bins of equal size. With Nb the number of bins, and Ns the number of multi-region samples. The number of input samples should be 2 or more (Ns > 1). NOTE: the input bin-values should be non-log values, and should be normalized, i.e. have an average of 1 per sample.

OUTPUT

1) LPAC_success    Binary variable with 1 if LPAC inference was succesfull and 0 if unsuccesfull
2) ploidies        Ns x 1 size vector containing the tumor ploidies
3) purities        Ns x 1 size vector containing the sample purities
4) CNH             Ns x 1 size vector containing the copy number heterogeneity
5) best_sample     Ns x 1 size binary vector in which the best sample used for multi-region inference is identified.
