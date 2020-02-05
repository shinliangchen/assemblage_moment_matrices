function constr = AMM_proj_SDP_constraints(pax,pby,pabxy,gamma_SDP,u,uni_mono)
%%%%%%
%
% AMM_proj_SDP_constraints.m generates a set of SDP constraints required
% for the AMM (in the projector form)
%
%  pax,pby,pabxy : the SDP variables for probabilities in a Bell scenario
%       gamma_SDP: the SDP variables for AMM
%               u: the unknown variables in AMM
%        uni_mono: the unique elements associated with the monomials of AMM
%
% author: Shin-Liang Chen
%%%%%%
[na, nx] = size(pax);
[nb, ny] = size(pby);

constr = [];
    
for y = 1:ny
    for b = 1:(nb-1)
        str_Bby = string(strcat('B_',num2str(b),'|',num2str(y)));
        
        for x = 1:nx
            for a = 1:na
                constr = [constr, pabxy{a,b,x,y} == u{uni_mono==str_Bby,a,x}];
                constr = [constr, pabxy{a,b,x,y} >= 0];
            end
        end
    end
end

for x = 1:nx
    for a = 1:na
        for y = 1:ny
            constr = [constr, pabxy{a,nb,x,y} >= 0];
        end
    end
end

for x = 1:nx
    sum_a_gamma_SDP{x} = 0;
    for a = 1:na
        constr = [constr, u{1,a,x} == pax{a,x}]; % elements corresponding to string('Id') is pax{a,x}
        sum_a_gamma_SDP{x} = sum_a_gamma_SDP{x} + gamma_SDP{a,x};
    end
end

for x = 1:nx
    for a = 1:na
        constr = [constr, gamma_SDP{a,x}>=0, (1/2).*(gamma_SDP{a,x}'+gamma_SDP{a,x})>=0];
        constr = [constr, pax{a,x} >= 0];
    end
end

for x = 2:nx
    constr = [constr, sum_a_gamma_SDP{x} == sum_a_gamma_SDP{1}];
end

end