clear all;
tic;

nx = 4;
ny = 3;
na = 2;
nb = 2;
level = 1;

%S = [string('Id') string('B_1|1') string('B_1|2') string('B_1|3')];
S = AMM_proj_gen_xlevel_seq(ny,nb,level);

[gamma_SDP,u,uni_mono] = AMM_proj_SDP_variables_sequence(nx,na,S);
[pax,pby,pabxy] = SDP_Variables_For_Pabxy(nx, ny, na, nb);

for x = 1:nx
    for y = 1:ny
        E{x,y} = pabxy{1,1,x,y} + pabxy{2,2,x,y} - pabxy{1,2,x,y} - pabxy{2,1,x,y};
        
    end
end


Bell_inequality = E{1,1} + E{1,2} + E{1,3} + E{2,1} - E{2,2} - E{2,3} + ...
     -E{3,1} + E{3,2} - E{3,3} - E{4,1} - E{4,2} + E{4,3};

constr = AMM_proj_SDP_constraints(pax,pby,pabxy,gamma_SDP,u,uni_mono);

sol = solvesdp(constr , -Bell_inequality);
sol
value_SE = double(Bell_inequality)

