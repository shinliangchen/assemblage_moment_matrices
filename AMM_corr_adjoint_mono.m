function adj_result = AMM_corr_adjoint_mono(monomial)
%%%%%%
% AMM_corr_adjoint_mono.m reverses the order of the product of observables
% E.g. transform B_2*B_3 to B_3*B_2.
%
% monomial: the correlator monomial in the string form. E.g., string('B_2*B_1').
%
% Note: We strongly recommend one to choose observables from the set {B_1,
% B_2,...,B_9}
% 
% author: Shin-Liang Chen
%%%%%%

if monomial == string('Id')
    
    adj_result = string('Id');
    
else
    
    monomial_sep = strsplit(monomial,'*');
    
    if length(monomial_sep)==1
        adj_result = monomial_sep;
    else
        monomial_sep_dag = flip(monomial_sep);
        monomial_sep_dag = strjoin(monomial_sep_dag,'*');
        adj_result = strjoin(string(monomial_sep_dag),'*');
    end
    
end

end