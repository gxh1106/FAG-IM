function [Index_pattern, Sym_pattern, Vector_pattern] = gen_patterns(Na, N, const)
Np = N/Na;          % the number of ports in each group
M = length(const);
%% Generate Tx Patterns
% FA index pattern
Index_pattern = zeros(Na, Np^Na);
for idx0 = 1:Na
    idx1 = 1;
    idx2 = 1;
    while idx1 <= Np^Na
        for n = 1:Np^(idx0-1)
            Index_pattern(idx0,idx1) =  idx2;
            idx1 = idx1+1;
        end
        if mod(idx2+1,Np) == 0
            idx2 = Np;
        else
            idx2 = mod(idx2+1,Np);
        end
    end
end
% II = 1:Np;

% constellation pattern
Sym_pattern = zeros(Na, M^Na);
for idx0 = 1:Na
    idx1 = 1;
    idx2 = 1;
    while idx1 <= M^Na
        for n = 1:M^(idx0-1)
            Sym_pattern(idx0,idx1) =  const(idx2);
            idx1 = idx1+1;
        end
        if mod(idx2+1,M) == 0
            idx2 = M;
        else
            idx2 = mod(idx2+1,M);
        end
    end
end

% FA index with constellation symbols pattern
Vector_pattern = zeros(N, Np^Na * M^Na);
for ii = 1:Np^Na
    Vector_pattern( Index_pattern(:,ii) + (0:Na-1).' .* Np, (ii-1)*M^Na + 1 : ii*M^Na) = Sym_pattern;
end

end