function string_moment_matrix = AMM_corr_string_real(Seq)
%%%%%%
%
% AMM_proj_string.m generates the representation of the assemblage moment
% matrices associated with the sequence S, under the assumption that the
% AMM are real.
%
%              Seq: a sequence of Bob's projector
%           Seqdag: the complex conjugate of S (use proj_adjoint_poly.m to generate)
%  real_or_complex: 'real' for real AMM ; 'compplex' for complex AMM
%
% author: Shin-Liang Chen
%%%%%%
for ii = 1:length(Seq)
    Seqdag(ii) = AMM_corr_adjoint_mono(Seq(ii));
end

string_moment_matrix = string(zeros(length(Seq)));

for i = 1:length(Seqdag)
    for j = 1:length(Seq)
        Sdag_i = strsplit(Seqdag(i),'*');
        S_j = strsplit(Seq(j),'*');
        Si_and_Sj = [Sdag_i S_j];
        Si_and_Sj_B = Si_and_Sj(~cellfun('isempty', strfind(Si_and_Sj,'B')));
        % this is to only keep moments of B
        
        if isempty(Si_and_Sj_B) % length of gamma_B = 0
            Si_and_Sj_B_short = string('Id');
        elseif length(Si_and_Sj_B) == 1
            Si_and_Sj_B_short = Si_and_Sj_B;
            % gamma_B_short = Bk, e.g., B3
        elseif length(Si_and_Sj_B) >= 1
            Si_and_Sj_B_compu = Si_and_Sj_B;
            for i_B = 2:length(Si_and_Sj_B)
                if Si_and_Sj_B_compu(i_B) == Si_and_Sj_B_compu(i_B-1)
                    Si_and_Sj_B_compu(i_B-1) = string('Id');
                    Si_and_Sj_B_compu(i_B) = string('Id');
                    % the two are combined to two identities if they are the same,
                    % e.g. B_1*B_1 = Id*Id
                    Si_and_Sj_B_compu = [Si_and_Sj_B_compu(cellfun('isempty', strfind(Si_and_Sj_B_compu,'B')))...
                        Si_and_Sj_B_compu(~cellfun('isempty', strfind(Si_and_Sj_B_compu,'B')))];
                    % this is to move Id to the most-left part, e.g., change [B_2 Id Id B_1] to [Id Id B_2 B_1]
                end
            end
            Si_and_Sj_B_short = Si_and_Sj_B_compu(cellfun('isempty', strfind(Si_and_Sj_B_compu,'Id'))); % remove Id
            
        else
            error('something is wrong')
        end
        
        Si_and_Sj_B_short_flip = flip(Si_and_Sj_B_short);
        
        if isempty(Si_and_Sj_B_short) == 1
            string_moment_matrix(i,j) = string('Id');
        elseif ~or(any(any(strjoin(Si_and_Sj_B_short,'*') == string_moment_matrix)) ,...
                any(any(strjoin(Si_and_Sj_B_short_flip,'*') == string_moment_matrix)))
            string_moment_matrix(i,j) = strjoin(Si_and_Sj_B_short,'*');
        elseif any(any(strjoin(Si_and_Sj_B_short,'*') == string_moment_matrix))
            string_moment_matrix(i,j) = strjoin(Si_and_Sj_B_short,'*');
        elseif any(any(strjoin(Si_and_Sj_B_short_flip,'*') == string_moment_matrix))
            string_moment_matrix(i,j) = strjoin(Si_and_Sj_B_short_flip,'*');
        else
            error('something is wrong')
        end
        
    end
end

end