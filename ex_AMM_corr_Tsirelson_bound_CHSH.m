clear all;
tic;

nx = 2;
ny = 2;
na = 2;
nb = 2;
level = 1;

%S = [string('Id') string('B_1') string('B_2')];
S = AMM_corr_gen_xlevel_seq(ny,1);

[gamma_SDP,u,uni_mono] = AMM_corr_SDP_variables_sequence(nx,na,S);
[pax,pby,pabxy] = SDP_Variables_For_Pabxy(nx, ny, na, nb);
[constr,E] = AMM_corr_SDP_constraints(pax,pby,pabxy,gamma_SDP,u,uni_mono);

Bell_inequality = E{1,1}+E{1,2}+E{2,1}-E{2,2};

sol = solvesdp(constr , -Bell_inequality);
sol
value_SE = double(Bell_inequality)
toc
