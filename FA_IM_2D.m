clear all; close all;

rng(0)

N_frame = 10;         % the number of test frames 100000 50000
SNR = 4 : 2 : 22;

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

M = 2;              % modulation order
% if M == 2
%     const = [1, -1];
% else
%     const = QAM_init(M);
% end
% Dem = comm.GeneralQAMDemodulator(const);
const = qammod(0:M-1, M, 'gray', 'UnitAveragePower',true);

%% Caculate the Spectral Efficiency

% NoGroup
A = nchoosek(N,G);                 % possible FA position patterns available
m_index = floor(log2(A));        % FA grids index bits
m_mod = log2(M);                    % constellation bits
bpcu = m_index + G*m_mod;    % bits per channel use (SE)

%% Generate Channel
N_chan = 100;          % paper: num_realization = 10000;

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

result_folder=['result/VS/'];
if ~exist(result_folder,'dir')
    mkdir(result_folder)
end
result_file = [result_folder,'FA_IM_bpcu_',num2str(bpcu),'_W',num2str(W1),'x',num2str(W2),'_Nr_',num2str(Nr),'.mat'];
config.N = NN;
config.W = W;
config.G = GG;
config.Nr = Nr;
config.M = M;

%% Generate Tx Patterns
% FA index pattern
K = 2^m_index;
II = nchoosek(1:N, G);
II = II([1:K/2, end-K/2+1:end], :);
% II = II(1:K, :);

% if Np>Na
%     II_a = nchoosek(1:Np, Na);
%     KK = size(II_a, 1);
%     II_aa = repmat(II_a, Na, 1) + kron((0:Na-1).' .* Np, ones(KK, Na));
%     II = II(1:K, :);
% end

% constellation pattern
sym_pattern = zeros(G, M^G);
for idx0 = 1:G
    idx1 = 1;
    idx2 = 1;
    while idx1 <= M^G
        for n = 1:M^(idx0-1)
            sym_pattern(idx0,idx1) =  const(idx2);
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
IIS_pattern = zeros(N, K * M^G);
for ii = 1:K
    IIS_pattern( II(ii,:), (ii-1)*M^G + 1 : ii*M^G) = sym_pattern;
end


%% Generate Tx Antenna Grids Index
N_ibit = m_index * N_frame;
ibit_src = randi([0 1], N_ibit, 1);   % bit sequence
ibit_mat = reshape(ibit_src, m_index, []).';
II_idx = bi2de(ibit_mat) + 1;
% I_idx = II(II_idx,:);

%% Generate Tx Symbols
N_sbit = G * m_mod * N_frame;
sbit_src = randi([0 1], N_sbit, 1);   % bit sequence
sbit_mat = reshape(sbit_src, m_mod, []).';
sym_idx = bi2de(sbit_mat)+1;
s = reshape(const(sym_idx), G, N_frame);  % Complex transmit symbols vector

x = zeros(N, N_frame);
for idx3 = 1 : N_frame
    x( II(II_idx(idx3), :),  idx3) = s(:, idx3);
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
    disp(['Iter -- ',num2str(idx4)])
    SNRdB = SNR(idx4);
    N0 = 10^(-SNRdB/10);
    % N0 = G * power * 10^(-SNRdB/10);  
    % Received Signal
    
    err_sbit = zeros(N_chan, 1);
    err_ibit = zeros(N_chan, 1);
    for idx5 = 1 : N_chan
    % parfor idx5 = 1 : N_reali
        wn = sqrt(N0/2) .* (randn(Nr, N_frame) + 1j*randn(Nr, N_frame));
        % Hx(:, (jj-1)*N_frame + 1 : jj*N_frame) = H(:,:, jj)*x;
        Heff = H(:,:, idx5);      % 后续加上信道排序！

        y = Heff*x + wn;
        decode_pattern = sqrt(power) .* Heff*IIS_pattern;

        sym_dec = zeros( G*N_frame, 1 );
        II_dec = zeros( N_frame, 1 );
        for idx6 = 1 : N_frame
            distance = sum( abs(repmat(y(:, idx6), 1, K * M^G) - decode_pattern) .^ 2 , 1);
            [~, iis_dec] = min(distance);

            sym_pattern_dec = mod(iis_dec - 1, M^G) + 1;
            sym_dec( (idx6-1)*G + 1 : idx6*G ) = sym_pattern(:, sym_pattern_dec);
            II_dec(idx6) = floor( (iis_dec - 1) / M^G );
        end
        % sym_idx_dec = step(Dem, sym_dec)+1;  % Complex transmit symbols
        sym_idx_dec = qamdemod(sym_dec, M, 'gray','UnitAveragePower', true)+1;
        sbit_mat_dec = de2bi(sym_idx_dec-1, m_mod);
        sbit_dec = reshape(sbit_mat_dec.', [], 1); 

        ibit_mat_dec = de2bi(II_dec, m_index);
        ibit_dec = reshape(ibit_mat_dec.', [], 1);

        err_sbit(idx5) = sum(xor(sbit_dec, sbit_src));
        err_ibit(idx5) = sum(xor(ibit_dec, ibit_src));
    end
    BER_ML(idx4) = (sum(err_sbit) + sum(err_ibit)) / (N_chan * (N_sbit+N_ibit));
    BER_ML_I(idx4) = sum(err_ibit) / (N_chan * N_ibit);
    BER_ML_S(idx4) = sum(err_sbit) / (N_chan * N_sbit);

    save(result_file,'SNR','BER_ML','BER_ML_I','BER_ML_S');
end
figure
semilogy(SNR, BER_ML, 'b--o'), hold on
% semilogy(SNR, BER_ML_I, 'b--x'), hold on
% semilogy(SNR, BER_ML_S, 'b--^'), hold on
% semilogy(SNR, BER_ZF, 'r-o'), hold on

grid on
xlabel('SNR [dB]')
ylabel('BER')
% legend('ML','Index Bits','Symbol Bits')

axis([0 30 1e-5 1e-1])

save(result_file,'SNR','BER_ML','BER_ML_I','BER_ML_S');
% delete(gcp)