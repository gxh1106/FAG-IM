clear all; close all;

rng(0)

group = 'BLKG'; % 'MCG', 'BLKG'

N_frame = 10;         % the number of test frames 100000 50000
SNR = 6 : 2 : 16;

Nr = 8;             % the number of receive antennas

W = [1.6 1.6];
NN = [4 4];
GG = [2 2];

W1 = W(1);
W2 = W(2);
N1 = NN(1);             % the number of ports along vertical direction
N2 = NN(2);             % the number of ports along horizontal direction
N = N1 * N2;            % the number of ports           
G1 = GG(1);             % the number of groups along vertical direction
G2 = GG(2);             % the number of groups along horizontal direction
G = G1*G2;              % the number of active ports = the number of groups

P = N/G;          % the number of ports in each group
P1 = N1/G1;       % the number of ports in each group along vertical direction
P2 = N2/G2;       % the number of ports in each group along horizontal direction

M = 4;              % modulation order
const = qammod(0:M-1, M, 'gray', 'UnitAveragePower',true);

%% Caculate the Spectral Efficiency

% Group
m_index = floor(log2(P));      % FA grids index bits
m_mod = log2(M);                    % constellation bits
bpcu = G* (m_index + m_mod);   % bits per channel use (SE)

%% Generate Channel
N_chan = 100;          % paper: num_realization = 10000;

%% index first along NN1 direction then NN2 (NN2 times NN1)
N = NN(1)*NN(2);
% Antennas parameters
d = W./(NN-1);
y_pos = (0:NN(2)-1)*d(2);
x_pos = (0:NN(1)-1)*d(1);
[Xpos, Ypos] = ndgrid(x_pos, y_pos);
xy_pos = [Xpos(:) Ypos(:)];
Sigma = SigmaIso3D(xy_pos);
[V,Lambda] = eig(Sigma);
% Sort eigenvalues in descending order
[lambda, index] = sort(diag(Lambda),'descend');
% Only the eigenvalues larger than a small tolerance are stored
% lambda = lambda(lambda>1e-5);
% Re-arranges the corresponding eigenvectors
V = V(:,index(1:length(lambda)));
R = diag(sqrt(lambda))'*V';

H = zeros(Nr, N, N_chan);
for idx_reali = 1:N_chan
    Hc = sqrt(1/2) .* (randn(Nr, N) + 1j*randn(Nr, N));    
    H(:,:, idx_reali) = sqrt(1/G) .* Hc * R;
end

% H = FAS_channel_new(G1, G2, P1, P2, W1, W2, Nr, N_chan);

% channel_folder='channel_file';
% channel_file = [channel_folder,'/W1_',num2str(W1),'_W2_',num2str(W2),'_N_',num2str(N),'_',num2str(N1),'x',num2str(N2),'_Na_',num2str(G),'_',num2str(G1),'x',num2str(G2),'_Nr_',num2str(Nr),'.mat'];
% H = load(channel_file).H;
% N_chan = size(H, 3); 

% Sigma2 = 1/G;
% H = sqrt(Sigma2) .* H;

%% grouping
if strcmp(group, 'BLKG')    
    index_H = zeros(N, N_chan);
    for idx_chan = 1:N_chan  
        H_sort_idx = zeros(P, G);
        for g = 1 : G
            g1 = mod(g-1, G1) + 1;
            g2 = ceil(g/G1);
            [p1, p2] = ndgrid((g2-1)*P2*N1 + (g1-1)*P1 + (1 : P1), (0:P2-1).*N1);
            H_sort_idx(:, g) = p1(:) +p2(:);
        end
        index_H(:,idx_chan) = H_sort_idx(:);
        H(:,:,idx_chan) = H(:, index_H(:,idx_chan), idx_chan);
    end
elseif strcmp(group, 'MCG')   
    index_H = zeros(N, N_chan);
    for idx_chan = 1:N_chan  
        H_tmp = H(:,:,idx_chan);
        port_table = 1:N;

        H_sort_idx = zeros(P, G);
        for g = 1 : G-1
            metric_dot = zeros(length(port_table) - 1, 1);
            for idx_dot = 1 : length(port_table) - 1
                metric_dot(idx_dot) = abs(dot( H_tmp(:, port_table(1)), H_tmp(:, port_table(idx_dot+1)) ));
            end
            [~, loc_dot] = maxk(metric_dot, P-1);
            port_idx = [1;loc_dot+1];
            H_sort_idx(:, g) = port_table(port_idx);
            port_table = port_table(setdiff(1 : length(port_table), port_idx));
        end
        H_sort_idx(:, G) = port_table;
        index_H(:,idx_chan) = H_sort_idx(:);
        H(:,:,idx_chan) = H(:, index_H(:,idx_chan), idx_chan);
    end
end


result_folder=['result/VS/'];
if ~exist(result_folder,'dir')
    mkdir(result_folder)
end
result_file = [result_folder,'FAG_IM_bpcu_',num2str(bpcu),'_W',num2str(W1),'x',num2str(W2),'_Nr_',num2str(Nr),'.mat'];
config.N = NN;
config.W = W;
config.G = GG;
config.Nr = Nr;
config.M = M;

[Index_pattern, Sym_pattern, Vector_pattern] = gen_patterns(G, N, const);

%% Generate Tx Antenna Grids Index
N_ibit = G * m_index * N_frame;
ibit_src = randi([0 1], N_ibit, 1);   % bit sequence
ibit_mat = reshape(ibit_src, m_index, []).';
II_idx = bi2de(ibit_mat) + 1;
II = reshape(II_idx, G, N_frame) + repmat((0:G-1).' .* P, 1, N_frame);  
% I_idx = II(II_idx,:);

%% Generate Tx Symbols
N_sbit = G * m_mod * N_frame;
sbit_src = randi([0 1], N_sbit, 1);   % bit sequence
sbit_mat = reshape(sbit_src, m_mod, []).';
sym_idx = bi2de(sbit_mat)+1;
s = reshape(const(sym_idx), G, N_frame);  % Complex transmit symbols vector

x = zeros(N, N_frame);
for idx3 = 1 : N_frame
    x( II(:, idx3),  idx3) = s(:, idx3);
end

power = 1;
x = sqrt(power) .* x;

% num_core = 90;
% parpool(num_core)

BER_ML = zeros(1, length(SNR));
BER_ML_I = zeros(1, length(SNR));
BER_ML_S = zeros(1, length(SNR));
BER_ZF = zeros(1, length(SNR));
BER_MMSE = zeros(1, length(SNR));
for idx4 = 1:length(SNR)
    disp(['Iter -- ',num2str(idx4), ' of ', num2str(length(SNR))])
    SNRdB = SNR(idx4);
    N0 = 10^(-SNRdB/10);
    % N0 = Na * power * 10^(-SNRdB/10);  
    % Received Signal
    
    err_sbit = zeros(N_chan, 1);
    err_ibit = zeros(N_chan, 1);
    for idx5 = 1 : N_chan
    % parfor idx5 = 1 : N_chan
        wn = sqrt(N0/2) .* (randn(Nr, N_frame) + 1j*randn(Nr, N_frame));
        Heff = H(:,:, idx5);

        y = Heff*x + wn;
        decode_pattern = sqrt(power) .* Heff*Vector_pattern;

        sym_dec = zeros( G*N_frame, 1 );
        II_dec = zeros( G*N_frame, 1 );
        for idx6 = 1 : N_frame
            distance = sum( abs(repmat(y(:, idx6), 1, P^G * M^G) - decode_pattern) .^ 2 , 1);
            [~, iis_dec] = min(distance);

            sym_pattern_dec = mod(iis_dec - 1, M^G) + 1;
            sym_dec( (idx6-1)*G + 1 : idx6*G ) = Sym_pattern(:, sym_pattern_dec);
            II_pattern_dec = ceil(iis_dec/(M^G));
            II_dec( (idx6-1)*G + 1 : idx6*G ) = Index_pattern(:, II_pattern_dec);
        end
        % sym_idx_dec = step(Dem, sym_dec)+1;  % Complex transmit symbols
        sym_idx_dec = qamdemod(sym_dec, M, 'gray','UnitAveragePower', true)+1;
        sbit_mat_dec = de2bi(sym_idx_dec-1, m_mod);
        sbit_dec = reshape(sbit_mat_dec.', [], 1);   

        ibit_mat_dec = de2bi(II_dec-1, m_index);
        ibit_dec = reshape(ibit_mat_dec.', [], 1);

        err_sbit(idx5) = sum(xor(sbit_dec, sbit_src));
        err_ibit(idx5) = sum(xor(ibit_dec, ibit_src));
    end
    BER_ML(idx4) = (sum(err_sbit) + sum(err_ibit)) / (N_chan * (N_sbit + N_ibit));
    BER_ML_I(idx4) = sum(err_ibit) / (N_chan * N_ibit);
    BER_ML_S(idx4) = sum(err_sbit) / (N_chan * N_sbit);

    save(result_file,'config','SNR','BER_ML','BER_ML_I','BER_ML_S');
end
figure
semilogy(SNR, BER_ML, 'r-o'), hold on
% semilogy(SNR, BER_ML_I, 'r-x'), hold on
% semilogy(SNR, BER_ML_S, 'r-^'), hold on
% semilogy(SNR, BER_ZF, 'r-o'), hold on

grid on
xlabel('SNR [dB]')
ylabel('BER')
% legend('ML','Index Bits','Symbol Bits')

axis([0 30 1e-5 1e-0])

save(result_file,'config','SNR','BER_ML','BER_ML_I','BER_ML_S');
% delete(gcp)