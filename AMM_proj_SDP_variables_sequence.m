function [gamma_SDP,u,uni_mono] = AMM_proj_SDP_variables_sequence(nx,na,Seq)
%%%%%%
%
% AMM_proj_SDP_variables_sequence.m generates a set of SDP variables (gamma_SDP
% and u) and the unique elements associated with the monomials of AMM
% (uni_mono)
%
%  gamma_SDP: the SDP variables for AMM
%          u: the unknown variables in AMM
%
%         nx: the number of Alice's measurement settings
%         na: the number of Alice's measurement outcomes
%        Seq: a sequence of Bob's projector
%
% author: Shin-Liang Chen
%%%%%%
gamma_str_real = AMM_proj_string_real(Seq);

uni_mono = unique(gamma_str_real);

for x = 1:nx
    for a = 1:na
        gamma_SDP{a,x} = zeros(length(gamma_str_real));
    end
end

uni_mono(uni_mono==string('0'))=[];
if any(uni_mono==string('Id'))
    uni_mono(uni_mono==string('Id')) = [];
    uni_mono = [string('Id'); uni_mono];
else
    disp('there is no identity operator in the AMM')
end

for a = 1:na
    for x = 1:nx
        for idx = 1:length(uni_mono)
            
            tfMatrix = (gamma_str_real==uni_mono(idx));
            u{idx,a,x} = sdpvar(1,1,'hermitian','real');
            gamma_SDP{a,x} = gamma_SDP{a,x} + u{idx,a,x}.*tfMatrix;
            
            if mod(idx,50)==0
                disp(strcat('generating SDP variables in loops of:', num2str(idx)))
            end
        end
    end
end

end