function adj_result = AMM_proj_adjoint_mono(monomial)
%%%%%%
% AMM_proj_adjoint_mono.m reverses the order of the product of projectors
% E.g. transform B_2|1*B_1|3 to B_1|3*B_2|1.
%
% monomial: the projector monomial in the string form. E.g., string('B_2|1*B_1|2').
%
% Note: We strongly recommend one to choose observables from the set
% {B_b|y} with b and y running from 1 to 9.
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