function [pax,pby,pabxy] = SDP_Variables_For_Pabxy(nx, ny, na, nb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SDP_Variables_For_Pabxy.m generates a set of SDP varialbes for 
% probability distributions P(a,b|x,y)
%
% nx: number of measurement settings for Alice
% ny: number of measurement settings for Bob
% na: number of measurement outcomes for Alice
% nb: number of measurement outcomes for Bob
%
% author: Shin-Liang Chen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for x = 1:nx
    
    sum_pax_except_last{x} = 0;
    
for a = 1:na-1
    
    pax{a,x} = sdpvar(1,1,'hermitian','real');
    sum_pax_except_last{x} = sum_pax_except_last{x} + pax{a,x};
    
end

    pax{na,x} = 1 - sum_pax_except_last{x};

end


for y = 1:ny
    
    sum_pby_except_last{y} = 0;
    
for b = 1:nb-1
    
    pby{b,y} = sdpvar(1,1,'hermitian','real');
    sum_pby_except_last{y} = sum_pby_except_last{y} + pby{b,y};
    
end

    pby{nb,y} = 1 - sum_pby_except_last{y};

end


for x = 1:nx
for y = 1:ny
    
    sum_ab_pabxy_except_last_ab{x,y} = 0;
    
for a = 1:na-1
    
    
for b = 1:nb-1
    
    pabxy{a,b,x,y} = sdpvar(1,1,'hermitian','real');
    sum_ab_pabxy_except_last_ab{x,y} = sum_ab_pabxy_except_last_ab{x,y} + pabxy{a,b,x,y};
    
end


end

end
end


for x = 1:nx
for y = 1:ny
    
    sum_a_pabxy_fix_b_in_last{x,y} = 0;
    
for a = 1:na-1
    
    sum_b_pabxy_except_last{a,x,y} = 0;
    
for b = 1:nb-1
    
    sum_b_pabxy_except_last{a,x,y} = sum_b_pabxy_except_last{a,x,y} + pabxy{a,b,x,y};
    
end

    pabxy{a,nb,x,y} = pax{a,x} - sum_b_pabxy_except_last{a,x,y};
    sum_a_pabxy_fix_b_in_last{x,y} = sum_a_pabxy_fix_b_in_last{x,y} + pabxy{a,nb,x,y};

end

end
end


for x = 1:nx
for y = 1:ny
    
    sum_b_pabxy_fix_a_in_last{x,y} = 0;
    
for b = 1:nb-1
    
    sum_a_pabxy_except_last{b,x,y} = 0;
    
for a = 1:na-1
    
    sum_a_pabxy_except_last{b,x,y} = sum_a_pabxy_except_last{b,x,y} + pabxy{a,b,x,y};
    
end

    pabxy{na,b,x,y} = pby{b,y} - sum_a_pabxy_except_last{b,x,y};
    sum_b_pabxy_fix_a_in_last{x,y} = sum_b_pabxy_fix_a_in_last{x,y} + pabxy{na,b,x,y};

end
end
end


for x = 1:nx
for y = 1:ny
    
    pabxy{na,nb,x,y} = 1 - sum_ab_pabxy_except_last_ab{x,y} - sum_a_pabxy_fix_b_in_last{x,y} - sum_b_pabxy_fix_a_in_last{x,y};

end
end


end