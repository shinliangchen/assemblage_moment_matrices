function seq_all = AMM_corr_gen_xlevel_seq(ny,level)
%%%%%%
%
% AMM_corr_gen_xlevel_seq.m generates a sequence of Bob's correlator with
% binary measurement outcomes
%
%      ny: the number of Bob's measurement settings
%   level: the level of moment relaxation
%
% author: Shin-Liang Chen
%%%%%%

if any(level == [1 2 3 4]) == 0
    error('for practical usage, the level should be set as a positive integer lower or equal than 4')
end

seq_S_B = [];
for y = 1:ny
    seq_S_B = [seq_S_B string(strcat('B_',num2str(y)))];
end

if level == 1
    seq_all = [string('Id') seq_S_B];
else
    seq_all = [string('Id') seq_S_B];
    seq_S_xB{1} = seq_S_B;
    
    for i_level = 2:level
        Sx = seq_S_xB{i_level-1};
        S1 = seq_S_xB{1};
        
        idx = 1;
        for i = 1:length(Sx)
            for j = 1:length(S1)
                Sx_i = strsplit(Sx(i),'*');
                S1_j = strsplit(S1(j),'*');
                Sxi_and_Sxj = [Sx_i S1_j];
                
                gamma_B = Sxi_and_Sxj;
                
                gamma_B_proj = gamma_B;
                for i_B = 2:length(gamma_B)
                    by1 = strsplit(gamma_B_proj(i_B-1),'_');
                    y1 = str2num(char(by1(2)));
                    by2 = strsplit(gamma_B_proj(i_B),'_');
                    y2 = str2num(char(by2(2)));
                    
                    if gamma_B_proj(i_B) == gamma_B_proj(i_B-1)
                        gamma_B_proj(i_B-1) = string('Id');
                        gamma_B_proj(i_B) = string('Id');
                        % the two are combined to two identities if they are the same,
                        % e.g. B_1*B_1 = Id*Id
                        gamma_B_proj = [gamma_B_proj(cellfun('isempty', strfind(gamma_B_proj,'B')))...
                            gamma_B_proj(~cellfun('isempty', strfind(gamma_B_proj,'B')))];
                        % this is to move Id to the most-left part, e.g., change [B_2 Id Id B_3] to [Id Id B_2 B_3]
                    end
                end
                
                gamma_B_short = gamma_B_proj(cellfun('isempty', strfind(gamma_B_proj,'Id')));
                
                if length(gamma_B_short)==i_level
                    seq_S_xB{i_level}(idx) = strjoin(gamma_B_short,'*');
                    idx = idx + 1;
                end
                
            end
        end
        seq_all = [seq_all, seq_S_xB{i_level}];
    end
    
end

end