% Here we numerically check, for the ideal case, the relations of the
% equality before Eq.(B12) and Eq.(B12). Namely,
% 1) the relation between -i*tr(rho_{a|x}*B_5*B_6*B_7) and q
% 2) the relation between -i*tr(rho_{a|x}*B_k*B_l) and q for
%    (x,k,l)=(1,6,7),(2,7,5),(3,5,6)

clear all;

ny = 4; % number of measurement settings
nx = 3;
% nb = 2; % number of outcomes
na = 2;
% n = 2;

pauliZ = [1 0;0 -1];
pauliX = [0 1;1 0];
pauliY = [0 -1i;1i 0];

A{1} = pauliZ;
A{2} = pauliX;
A{3} = pauliY;



for x = 1:nx
    EA{1,x} = (1/2).*(eye(2)+A{x}); % +1 for the 1st outcome
    EA{2,x} = (1/2).*(eye(2)-A{x}); % -1 for the 2nd outcome
end

B{1} = (1/sqrt(3)).*(pauliZ+pauliX-pauliY);
B{2} = (1/sqrt(3)).*(pauliZ-pauliX+pauliY);
B{3} = (1/sqrt(3)).*(-pauliZ+pauliX+pauliY);
B{4} = (1/sqrt(3)).*(-pauliZ-pauliX-pauliY);


for y = 1:ny
    EB{1,y} = (1/2).*(eye(2)+B{y});
    EB{2,y} = (1/2).*(eye(2)-B{y});
end

rhoAB = (1/2).*[1 0 0 1;0 0 0 0;0 0 0 0;1 0 0 1];

state0 = (1/2).*(eye(2)+pauliZ);
state1 = (1/2).*(eye(2)-pauliZ);

n_p = 100;
if n_p > 10000
    error('n_p may be too large')
end

qqq = linspace(0,1,n_p);
tic;
for i_p = 1:n_p

qq = qqq(i_p);
for x = 1:nx
    for a = 1:na
        sigma_ax(:,:,a,x) = TrX(kron(EA{a,x},eye(2))*rhoAB,1,[2 2]);
        sigma_ax_ref(:,:,a,x) = qq.*kron(sigma_ax(:,:,a,x),state0) + (1-qq).*kron(transpose(sigma_ax(:,:,a,x)),state1);
        % this is the controlled mixture
    end
end


B5 = kron(pauliZ,eye(2));
B6 = kron(pauliX,eye(2));
B7 = kron(pauliY,pauliZ);


%%%%% check the quality before Eq.(B12) (start) %%%%%%
tr_B_rhoax_56711(i_p) = -1i.*trace(sigma_ax_ref(:,:,1,1)*(B5*B6*B7));
tr_B_rhoax_56721(i_p) = -1i.*trace(sigma_ax_ref(:,:,2,1)*(B5*B6*B7));
tr_B_rhoax_56712(i_p) = -1i.*trace(sigma_ax_ref(:,:,1,2)*(B5*B6*B7));
tr_B_rhoax_56722(i_p) = -1i.*trace(sigma_ax_ref(:,:,2,2)*(B5*B6*B7));
tr_B_rhoax_56713(i_p) = -1i.*trace(sigma_ax_ref(:,:,1,3)*(B5*B6*B7));
tr_B_rhoax_56723(i_p) = -1i.*trace(sigma_ax_ref(:,:,2,3)*(B5*B6*B7));

diff_Eq_B12_a(1,i_p) = tr_B_rhoax_56711(i_p) - (qq-1/2);
% check if tr_B_rhoax_56711 is equal to qq-1/2 (similar for the following 5
% terms)
diff_Eq_B12_a(2,i_p) = tr_B_rhoax_56721(i_p) - (qq-1/2);
diff_Eq_B12_a(3,i_p) = tr_B_rhoax_56712(i_p) - (qq-1/2);
diff_Eq_B12_a(4,i_p) = tr_B_rhoax_56722(i_p) - (qq-1/2);
diff_Eq_B12_a(5,i_p) = tr_B_rhoax_56713(i_p) - (qq-1/2);
diff_Eq_B12_a(6,i_p) = tr_B_rhoax_56723(i_p) - (qq-1/2);
%%%%% check the quality before Eq.(B12) (end) %%%%%%


%%%%% check Eq.(B12) (start) %%%%%%
for x = 1:nx
    for a = 1:na
        tr_B_rhoax_67(a,x,i_p) = -1i.*trace(sigma_ax_ref(:,:,a,x)*(B6*B7));
        tr_B_rhoax_75(a,x,i_p) = -1i.*trace(sigma_ax_ref(:,:,a,x)*(B7*B5));
        tr_B_rhoax_56(a,x,i_p) = -1i.*trace(sigma_ax_ref(:,:,a,x)*(B5*B6));
    end
end

for a = 1:na
    diff_Eq_B12_b_x1_B67(a,i_p) = tr_B_rhoax_67(a,1,i_p) - ((-1)^(a-1))*(qq-1/2);
    diff_Eq_B12_b_x2_B75(a,i_p) = tr_B_rhoax_75(a,2,i_p) - ((-1)^(a-1))*(qq-1/2);
    diff_Eq_B12_b_x3_B56(a,i_p) = tr_B_rhoax_56(a,3,i_p) - ((-1)^a)*(qq-1/2);
end
%%%%% check Eq.(B12) (end) %%%%%%


end
toc;

% The following results are both zero if the equality before Eq.(B12) holds
max(max(diff_Eq_B12_a))
min(min(diff_Eq_B12_a))



% The following results are all zero if Eq.(B12) holds 
max(max(diff_Eq_B12_b_x1_B67)) % (x,k,l) = (1,6,7)
min(min(diff_Eq_B12_b_x1_B67)) % (x,k,l) = (1,6,7)

max(max(diff_Eq_B12_b_x2_B75)) % (x,k,l) = (2,7,5)
min(min(diff_Eq_B12_b_x2_B75)) % (x,k,l) = (2,7,5)

max(max(diff_Eq_B12_b_x3_B56)) % (x,k,l) = (3,5,6)
min(min(diff_Eq_B12_b_x3_B56)) % (x,k,l) = (3,5,6)


