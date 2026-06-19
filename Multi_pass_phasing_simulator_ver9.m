%% Bellhop データの読み込み
[Arr, Pos] = read_arrivals_asc('JamTank_Arr.arr');
BER_matrix  = zeros(16,4);
OSNR_matrix = zeros(16,4);

PTR_BER_vector     = zeros(4,1);
PTR_OSNR_vector    = zeros(4,1);
PTRDFE_BER_vector  = zeros(4,1);
PTRDFE_OSNR_vector = zeros(4,1);

idx_r = 1;

%% DFEパラメータ（RLS版）
% q関数のサイドローブ構造から設定
% lag=-1 に 0.7 の強いサイドローブ → Nf はそれをカバーする
% lag=-25 付近にも 0.3 程度 → Nf >= 30 が必要
Nf      = 35;    % フィードフォワードタップ数
Nb      = 10;    % フィードバックタップ数
lambda  = 0.999; % RLS忘却係数（1に近いほど過去データを重視）
delta   = 0.01;  % RLS正則化（逆行列の初期値スケール）
Ntrain  = 400;   % トレーニングシンボル数（RLSはLMSより少なくて済む）

for idx_sd = 1:4

    %% エリア①：送信信号の生成
    Rs = 2000; Sps = 16; Fs = Rs * Sps; numSymbols = 2000;

    bits       = randi([0 1], 2 * numSymbols, 1);
    inPhase    = 2*bits(1:2:end) - 1;
    quadrature = 2*bits(2:2:end) - 1;
    symbols    = inPhase + 1i*quadrature;
    sigPower   = mean(abs(symbols).^2);

    txSignal = zeros(numSymbols * Sps, 1);
    for i = 1:numSymbols
        txSignal((i-1)*Sps + 1 : i*Sps) = symbols(i);
    end

    %% 全受信機共通の最小遅延
    global_min_delay = inf;
    for idx_rd_pre = 1:16
        d_pre     = Arr(idx_r, idx_rd_pre, idx_sd).delay;
        d_pre     = real(d_pre(real(d_pre) > 0));
        if ~isempty(d_pre)
            global_min_delay = min(global_min_delay, min(d_pre));
        end
    end

    ptr_buf_len      = length(txSignal) * 4;
    total_ptr_signal = zeros(ptr_buf_len, 1);
    total_autocorr   = zeros(ptr_buf_len, 1);

    for idx_rd = 1:16

        %% エリア②：チャネル適用・ノイズ付加
        amplitudes = Arr(idx_r, idx_rd, idx_sd).A;
        delays     = real(Arr(idx_r, idx_rd, idx_sd).delay);

        valid_idx  = (delays > 0);
        amplitudes = amplitudes(valid_idx);
        delays     = delays(valid_idx);
        delays     = delays - global_min_delay;

        valid_delay = (delays >= 0);
        amplitudes  = amplitudes(valid_delay);
        delays      = delays(valid_delay);
        if isempty(delays), continue; end

        maxDelaySamples = round(max(delays) * Fs);
        rxSignalBase    = zeros(length(txSignal) + maxDelaySamples, 1);
        for p = 1:length(delays)
            delaySamples = round(delays(p) * Fs);
            idx          = delaySamples + (1:length(txSignal));
            rxSignalBase(idx) = rxSignalBase(idx) + amplitudes(p) * txSignal;
        end

        SNR_dB      = 20;
        sigPower_ch = mean(abs(rxSignalBase).^2);
        noisePower  = sigPower_ch / (10^(SNR_dB / 10));
        noise       = sqrt(noisePower/2) * (randn(size(rxSignalBase)) + 1i*randn(size(rxSignalBase)));
        rxSignal    = rxSignalBase + noise;

        %% --- A. 通常の受信処理（素のBER & OSNR）---
        sampleOffset = round(Sps / 2);
        rxSymbols    = rxSignal(sampleOffset : Sps : sampleOffset + Sps*(numSymbols-1));

        alpha_single   = (symbols' * rxSymbols) / (symbols' * symbols);
        rxSymbols_comp = rxSymbols / alpha_single;

        errPower_single             = mean(abs(symbols - rxSymbols_comp).^2);
        OSNR_matrix(idx_rd, idx_sd) = 10 * log10(sigPower / errPower_single);

        dec_bits          = zeros(2 * numSymbols, 1);
        dec_bits(1:2:end) = real(rxSymbols_comp) > 0;
        dec_bits(2:2:end) = imag(rxSymbols_comp) > 0;
        BER_matrix(idx_rd, idx_sd) = sum(bits ~= dec_bits) / length(bits);

        %% --- B. PTR処理 ---
        delay_samples    = round(delays * Fs) + 1;
        h                = zeros(max(delay_samples), 1);
        h(delay_samples) = amplitudes;

        h_ptr_filter = conj(flipud(h));
        z_m          = conv(rxSignal, h_ptr_filter);

        len_zm = min(length(z_m), ptr_buf_len);
        total_ptr_signal(1:len_zm) = total_ptr_signal(1:len_zm) + z_m(1:len_zm);

        h_ac   = conv(h, conj(flipud(h)));
        len_ac = min(length(h_ac), ptr_buf_len);
        total_autocorr(1:len_ac) = total_autocorr(1:len_ac) + h_ac(1:len_ac);

    end % ◀ 受信機ループ

    %% エリア③：PTR合成信号の復調
    [~, peak_idx] = max(abs(total_autocorr));
    ptr_offset    = peak_idx + round(Sps/2) - 1;

    ptrSymbols = total_ptr_signal(ptr_offset : Sps : ptr_offset + Sps*(numSymbols-1));

    alpha_ptr       = (symbols' * ptrSymbols) / (symbols' * symbols);
    ptrSymbols_comp = ptrSymbols / alpha_ptr;

    errPower_ptr            = mean(abs(symbols - ptrSymbols_comp).^2);
    PTR_OSNR_vector(idx_sd) = 10 * log10(sigPower / errPower_ptr);

    dec_bits_ptr          = zeros(2 * numSymbols, 1);
    dec_bits_ptr(1:2:end) = real(ptrSymbols_comp) > 0;
    dec_bits_ptr(2:2:end) = imag(ptrSymbols_comp) > 0;
    PTR_BER_vector(idx_sd) = sum(bits ~= dec_bits_ptr) / length(bits);

    %% --- C. PTR + RLS-DFE処理 ---
    % RLSはLMSより高速に収束（収束に必要な反復≈ Nf+Nb 回）
    dfe_input = ptrSymbols;   % alpha補正前の生PTRシンボル

    Ntap = Nf + Nb;
    w    = zeros(Ntap, 1);    % [FF係数; FB係数]
    w(1) = 1;                 % 初期値：最初のFFタップを1

    % RLS逆相関行列の初期化
    P    = (1/delta) * eye(Ntap);

    ff_buf = zeros(Nf, 1);
    fb_buf = zeros(Nb, 1);

    dfe_symbols = zeros(numSymbols, 1);

    for n = 1:numSymbols
        ff_buf = [dfe_input(n); ff_buf(1:end-1)];

        % 入力ベクトル（FF部 + FB部）
        u = [ff_buf; -fb_buf];   % FBは減算なのでマイナス

        % DFE出力
        y = w' * u;
        dfe_symbols(n) = y;

        % QPSK硬判定
        d = sign(real(y)) + 1i * sign(imag(y));

        % 参照シンボル（トレーニング or 判定）
        if n <= Ntrain
            s_ref  = symbols(n);
            fb_new = symbols(n);
        else
            s_ref  = d;
            fb_new = d;
        end

        % RLS更新
        Pu    = P * u;
        kappa = Pu / (lambda + u' * Pu);   % カルマンゲイン
        e     = s_ref - y;                  % 誤差信号
        w     = w + kappa * conj(e);        % 係数更新
        P     = (P - kappa * (Pu)') / lambda; % 逆相関行列更新

        fb_buf = [fb_new; fb_buf(1:end-1)];
    end

    % DFE出力のOSNR・BER計算
    alpha_dfe        = (symbols' * dfe_symbols) / (symbols' * symbols);
    dfe_symbols_comp = dfe_symbols / alpha_dfe;

    errPower_dfe               = mean(abs(symbols - dfe_symbols_comp).^2);
    PTRDFE_OSNR_vector(idx_sd) = 10 * log10(sigPower / errPower_dfe);

    dec_bits_dfe          = zeros(2 * numSymbols, 1);
    dec_bits_dfe(1:2:end) = real(dfe_symbols_comp) > 0;
    dec_bits_dfe(2:2:end) = imag(dfe_symbols_comp) > 0;
    PTRDFE_BER_vector(idx_sd) = sum(bits ~= dec_bits_dfe) / length(bits);

end % ◀ 音源ループ

%% 選択ダイバーシティ
[selected_BER,  best_rx_ber_idx]  = min(BER_matrix, [], 1);
[selected_OSNR, best_rx_osnr_idx] = max(OSNR_matrix, [], 1);

%% エリア④：最終結果の出力
clc;
disp('==========================================================================');
disp(' 比較：【選択ダイバーシティ】 vs 【PTR (16ch合成)】 vs 【PTR + RLS-DFE】');
disp('==========================================================================');
for idx_sd = 1:4
    fprintf('【音源 %d】\n', idx_sd);
    fprintf('  ➔ 選択ダイバーシティ (Best 1ch性能):\n');
    fprintf('     - 最良BER : %f (受信機 %d)\n',      selected_BER(idx_sd),  best_rx_ber_idx(idx_sd));
    fprintf('     - 最高OSNR: %.2f dB (受信機 %d)\n', selected_OSNR(idx_sd), best_rx_osnr_idx(idx_sd));
    fprintf('  ➔ PTR（16chまるごとアレイ合成):\n');
    fprintf('     - 合成BER : %f\n',    PTR_BER_vector(idx_sd));
    fprintf('     - 合成OSNR: %.2f dB\n', PTR_OSNR_vector(idx_sd));
    fprintf('  ➔ PTR + RLS-DFE (Nf=%d, Nb=%d, Ntrain=%d):\n', Nf, Nb, Ntrain);
    fprintf('     - 合成BER : %f\n',    PTRDFE_BER_vector(idx_sd));
    fprintf('     - 合成OSNR: %.2f dB\n\n', PTRDFE_OSNR_vector(idx_sd));
end
