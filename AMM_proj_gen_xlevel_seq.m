function seq_all = AMM_proj_gen_xlevel_seq(ny,nb,level)
%%%%%%
%
% AMM_proj_gen_xlevel_seq.m generates a sequence of Bob's projectors with the
% form of Chen, Budroni, Liang, and Chen [Phys. Rev. Lett. 116, 240401 (2016)].
%
%      ny: number of Bob's measurement settings
%      nb: number of Bob's measurement outcomes
%   level: the level of moment relaxation
% 
% author: Shin-Liang Chen
%%%%%%

if any(level == [1 2 3 4]) == 0
    error('for practical usage, the level should be set as a positive integer lower or equal than 4')
end

seq_S_B = [];
for y = 1:ny
    for b = 1:nb-1
        seq_S_B = [seq_S_B string(strcat('B_',num2str(b),'|',num2str(y)))];
    end
end

if level == 1
    seq_all = [string('Id') seq_S_B];
else
    seq_all = [string('Id') seq_S_B];
    seq_S_xB{1} = seq_S_B;
    
    for i_level = 2:level
        Sdag = seq_S_xB{i_level-1};
        S = seq_S_xB{1};
        
        idx = 1;
        for i = 1:length(Sdag)
            for j = 1:length(S)
                Sdag_i = strsplit(Sdag(i),'*');
                S_j = strsplit(S(j),'*');
                gamma = [Sdag_i S_j];
                
                gamma_B = gamma(~cellfun('isempty', strfind(gamma,'B')));
                % this is to only keep moments of B
                is_zero = false;
                
                if isempty(gamma_B) % length of gamma_B = 0
                    gamma_B_short = string('Id');
                elseif length(gamma_B) == 1
                    gamma_B_short = gamma_B;
                    % gamma_B_short = B_b|y, e.g., B_1|2
                elseif length(gamma_B) >= 1
                    gamma_B_proj = gamma_B;
                    for i_B = 2:length(gamma_B)
                        by1 = strsplit(gamma_B_proj(i_B-1),'_');
                        by1 = strsplit(by1(2),'|');
                        b1 = str2num(char(by1(1)));
                        y1 = str2num(char(by1(2)));
                        by2 = strsplit(gamma_B_proj(i_B),'_');
                        by2 = strsplit(by2(2),'|');
                        b2 = str2num(char(by2(1)));
                        y2 = str2num(char(by2(2)));
                        
                        if and(y2==y1,b2~=b1)
                            is_zero = true;
                            break % stop the local loop
                        elseif gamma_B_proj(i_B) == gamma_B_proj(i_B-1)
                            gamma_B_proj(i_B) = string('Id');
                            % the two are combined to one if they are the same,
                            % e.g. B_1|1*B_1|1 = B_1|1*Id
                            gamma_B_proj = [gamma_B_proj(cellfun('isempty', strfind(gamma_B_proj,'B')))...
                                gamma_B_proj(~cellfun('isempty', strfind(gamma_B_proj,'B')))];
                            % this is to move Id to the most-left part, e.g., change [B_2|1 Id Id B_1|3] to [Id Id B_2|1 B_1|3]
                        end
                    end
                    gamma_B_short = gamma_B_proj(cellfun('isempty', strfind(gamma_B_proj,'Id')));
                else
                    error('something is wrong')
                end
                
                if length(gamma_B_short)==i_level
                    
                    gamma_B_short(gamma_B_short == string('Id'))=[];
                    gamma_temp = gamma_B_short;
                    
                    if isempty(gamma_temp)==1
                        
                    elseif is_zero == 1
                        
                    else
                        seq_S_xB{i_level}(idx) = strjoin(gamma_temp,'*');
                        idx = idx + 1;
                    end
                    
                end
                
            end
        end
        seq_all = [seq_all, seq_S_xB{i_level}];
    end
    
end

end
