function string_moment_matrix = AMM_proj_string_complex(Seq)
%%%%%%
%
% AMM_proj_string.m generates the representation of the assemblage moment
% matrices associated with the sequence S
%
%              Seq: a sequence of Bob's projector
%  real_or_complex: 'real' for real AMM ; 'compplex' for complex AMM
%
% author: Shin-Liang Chen
%%%%%%
for ii = 1:length(Seq)
    Seqdag(ii) = AMM_proj_adjoint_mono(Seq(ii));
end

string_moment_matrix = string(zeros(length(Seq)));

for i = 1:length(Seqdag)
    for j = 1:length(Seq)
        Sdag_i = strsplit(Seqdag(i),'*');
        S_j = strsplit(Seq(j),'*');
        Si_and_Sj = [Sdag_i S_j];
        Si_and_Sj_B = Si_and_Sj(~cellfun('isempty', strfind(Si_and_Sj,'B')));
        % The above is make a the product of moments be an array of
        % monomials (only keeping Bob's moments).
        is_zero = false;
        
        if isempty(Si_and_Sj_B)
            Si_and_Sj_B_B_short = string('Id');
        elseif length(Si_and_Sj_B) == 1
            Si_and_Sj_B_B_short = Si_and_Sj_B;
        elseif length(Si_and_Sj_B) >= 1
            % the following is to compute the product of S_j and Sdag_i
            Si_and_Sj_B_compu = Si_and_Sj_B;
            for i_B = 2:length(Si_and_Sj_B)
                by1 = strsplit(Si_and_Sj_B_compu(i_B-1),'_');
                by1 = strsplit(by1(2),'|');
                b1 = str2num(char(by1(1)));
                y1 = str2num(char(by1(2)));
                by2 = strsplit(Si_and_Sj_B_compu(i_B),'_');
                by2 = strsplit(by2(2),'|');
                b2 = str2num(char(by2(1)));
                y2 = str2num(char(by2(2)));
                
                if and(y2==y1,b2~=b1)
                    is_zero = true;
                    break % stop the local loop
                elseif Si_and_Sj_B_compu(i_B) == Si_and_Sj_B_compu(i_B-1)
                    Si_and_Sj_B_compu(i_B) = string('Id');
                    % the two are combined to one if they are the same,
                    % e.g. B_1|1*B_1|1 = B_1|1*Id
                    Si_and_Sj_B_compu = [Si_and_Sj_B_compu(cellfun('isempty', strfind(Si_and_Sj_B_compu,'B')))...
                        Si_and_Sj_B_compu(~cellfun('isempty', strfind(Si_and_Sj_B_compu,'B')))];
                    % this is to move Id to the most-left part, e.g., change [B_2|1 Id Id B_1|3] to [Id Id B_2|1 B_1|3]
                end
            end
            Si_and_Sj_B_B_short = Si_and_Sj_B_compu(cellfun('isempty', strfind(Si_and_Sj_B_compu,'Id')));
        else
            error('something is wrong')
        end
        
        if is_zero == 1
            string_moment_matrix(i,j) = string('0');
        else
            string_moment_matrix(i,j) = strjoin(Si_and_Sj_B_B_short,'*');
        end
        
    end
end

end