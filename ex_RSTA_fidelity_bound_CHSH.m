clear all;
%tic;

%%%%%% set parameters (start) %%%%%%
Bell_violation = 2*sqrt(2);

nx = 2;
ny = 2;
na = 2;
nb = 2;
level = 2;

c = cos(pi/8);
s = sin(pi/8);

S = AMM_corr_gen_xlevel_seq(ny,level);
%S = [string('Id') string('B_1') string('B_2')];% string('B_2*B_1') string('B_1*B_2')];

[gamma_SDP,u,uni_mono] = AMM_corr_SDP_variables_sequence(nx,na,S);
[pax,pby,pabxy] = SDP_Variables_For_Pabxy(nx, ny, na, nb);
[constr,E] = AMM_corr_SDP_constraints(pax,pby,pabxy,gamma_SDP,u,uni_mono);

SE = E{1,1}+E{1,2}+E{2,1}-E{2,2};

str_By1 = string('B_2*B_1');
if any(uni_mono==AMM_corr_adjoint_mono(str_By1))
    str_By1 = AMM_corr_adjoint_mono(str_By1);
elseif any(uni_mono==str_By1)==0
    error('something is wrong')
end

fidelity_all = (1/4).*...
    ( pax{1,1} + (c^2-s^2).*u{uni_mono==string('B_1'),1,1} + 2*c*s.*u{uni_mono==string('B_2'),1,1} - c*s.*u{uni_mono==str_By1,1,1} - c*s.*u{uni_mono==str_By1,1,1} +...
      pax{2,1} + (s^2-c^2).*u{uni_mono==string('B_1'),2,1} - 2*c*s.*u{uni_mono==string('B_2'),2,1} + c*s.*u{uni_mono==str_By1,2,1} + c*s.*u{uni_mono==str_By1,2,1} +...
      pax{1,2} + (c^2-s^2).*u{uni_mono==string('B_1'),1,2} - 2*c*s.*u{uni_mono==string('B_2'),1,2} + c*s.*u{uni_mono==str_By1,1,2} + c*s.*u{uni_mono==str_By1,1,2} +...
      pax{2,2} + (s^2-c^2).*u{uni_mono==string('B_1'),2,2} + 2*c*s.*u{uni_mono==string('B_2'),2,2} - c*s.*u{uni_mono==str_By1,2,2} - c*s.*u{uni_mono==str_By1,2,2} );


constr = [constr, SE==Bell_violation];

sol = solvesdp(constr, fidelity_all);
sol
double(fidelity_all)

%toc