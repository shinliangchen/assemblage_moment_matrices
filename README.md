# assemblage_moment_matrices

Before running the codes, make sure you have packages for solving semi-definite programs (SDP).

This project is devided into two parts:

 ## 1. The assemblage moment matrices (AMM)
 This part contains the related MATLAB codes of the AMM, which was proposed by Chen et al. [Phys. Rev. Lett. 116, 240401 (2018)] (see https://arxiv.org/abs/1603.08532 for the arXiv preprint).
 
 First, download all the codes and add all the files to the working path.
 
 Files with filenames started with "ex_..." are examples.
 
 For instance, "ex_AMM_corr_Tsirelson_bound_CHSH.m" computes an upper bound on the Tsirelson bound for a given CHSH inequality violation, using the AMM approach. In the filename, "corr" indicates that the AMM are in the form of correlators.
 
 The other example is "ex_AMM_proj_Tsirelson_bound_elegant_BI.m", which computes an upper bound on the Tsirelson bound for a given violation of the elegant Bell inequality. In the filename, "proj" indicates that the AMM are in the form of projectors, i.e., the form that originally proposed in the paper https://arxiv.org/abs/1603.08532 .
 
 In the above two examples, one can easily change them into any Bell scenario by adjusting the following parameters:
 nx - the number of Alice's measurement settings
 ny - the number of Bob's measurement settings
 na - the number of Alice's measurement outcomes
 nb - the number of Bob's measurement outcomes
 Bell_inequality - the Bell inequality under consideration
 
 ## 2. Robust self-testing of quantum assemblages
 This part includes the related MATLAB codes that computes a lower bound on the fidelity between the underlying assemblage and the reference one for a given Bell inequality violation.
 Ref: https://arxiv.org/abs/20xx.xxxxx
 
 For instance, "ex_RSTA_fidelity_bound_CHSH.m" computes a lower bound on the fidelity between the underlying assemblage and the CHSH-type assemblage for a given CHSH inequality violation.
 
