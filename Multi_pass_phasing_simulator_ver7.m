%% Bellhop データの読み込み
[Arr, Pos] = read_arrivals_asc('JamTank_Arr.arr');
BER_matrix  = zeros(16,4);
OSNR_matrix = zeros(16,4);

PTR_BER_vector     = zeros(4,1);
PTR_OSNR_vector    = zeros(4,1);
PTRDFE_BER_vector  = zeros(4,1);
PTRDFE_OSNR_vector = zeros(4,1);

idx_r = 1;

%% DFEパラメータ
Nf     = 16;    % フィードフォワードタップ数
Nb     = 8;     % フィードバックタップ数
mu_dfe = 0.005; % LMSステップサイズ
Ntrain = 200;   % トレーニングシンボル数

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

    % PTR用バッファ（十分大きく確保）
    ptr_buf_len      = length(txSignal) * 4;
    total_ptr_signal = zeros(ptr_buf_len, 1);

    % ★修正①：全受信機の自己相関を合算してピーク位置を正確に求めるバッファ
    total_autocorr = zeros(ptr_buf_len, 1);

    for idx_rd = 1:16

        %% エリア②：チャネル適用・ノイズ付加
        amplitudes = Arr(idx_r, idx_rd, idx_sd).A;
        delays     = Arr(idx_r, idx_rd, idx_sd).delay;

        valid_idx  = (delays > 0);
        amplitudes = amplitudes(valid_idx);
        delays     = delays(valid_idx);
        delays     = delays - min(delays);

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
        noise       = sqrt(noisePower / 2) * (randn(size(rxSignalBase)) + 1i*randn(size(rxSignalBase)));
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

        %% --- B. PTR処理（各受信機の寄与を合算）---
        delay_samples = round(delays * Fs) + 1;
        h             = zeros(max(delay_samples), 1);
        h(delay_samples) = amplitudes;

        h_ptr_filter = conj(flipud(h));
        z_m          = conv(rxSignal, h_ptr_filter);

        len_zm = length(z_m);
        if len_zm > ptr_buf_len
            z_m      = z_m(1:ptr_buf_len);
            len_zm   = ptr_buf_len;
        end
        total_ptr_signal(1:len_zm) = total_ptr_signal(1:len_zm) + z_m;

        % ★修正①：各受信機の自己相関を合算
        h_ac   = conv(h, conj(flipud(h)));
        len_ac = length(h_ac);
        if len_ac > ptr_buf_len
            h_ac   = h_ac(1:ptr_buf_len);
            len_ac = ptr_buf_len;
        end
        total_autocorr(1:len_ac) = total_autocorr(1:len_ac) + h_ac;

    end % ◀ 受信機ループの終わり

    %% エリア③：PTR合成信号の復調とBER/OSNR計算

    % ★修正①：全受信機の自己相関合算からピーク位置を決定
    [~, peak_idx] = max(abs(total_autocorr));
    ptr_offset    = peak_idx + round(Sps / 2) - 1;

    % シンボルサンプリング（インデックス範囲チェック付き）
    last_idx = ptr_offset + Sps * (numSymbols - 1);
    if last_idx > length(total_ptr_signal)
        warning('PTR信号のバッファが不足しています。ptr_buf_lenを増やしてください。');
    end
    ptrSymbols = total_ptr_signal(ptr_offset : Sps : ptr_offset + Sps*(numSymbols-1));

    % PTRのOSNR・BER（alpha補正）
    alpha_ptr       = (symbols' * ptrSymbols) / (symbols' * symbols);
    ptrSymbols_comp = ptrSymbols / alpha_ptr;

    errPower_ptr            = mean(abs(symbols - ptrSymbols_comp).^2);
    PTR_OSNR_vector(idx_sd) = 10 * log10(sigPower / errPower_ptr);

    dec_bits_ptr          = zeros(2 * numSymbols, 1);
    dec_bits_ptr(1:2:end) = real(ptrSymbols_comp) > 0;
    dec_bits_ptr(2:2:end) = imag(ptrSymbols_comp) > 0;
    PTR_BER_vector(idx_sd) = sum(bits ~= dec_bits_ptr) / length(bits);

    %% --- C. PTR + DFE処理 ---
    % ★修正③：DFEの入力はalpha補正前のptrSymbols（生シンボル）
    dfe_input = ptrSymbols;

    w_ff    = zeros(Nf, 1);
    w_ff(1) = 1;
    w_fb    = zeros(Nb, 1);
    ff_buf  = zeros(Nf, 1);
    fb_buf  = zeros(Nb, 1);

    dfe_symbols = zeros(numSymbols, 1);

    for n = 1:numSymbols
        ff_buf = [dfe_input(n); ff_buf(1:end-1)];

        y = w_ff' * ff_buf - w_fb' * fb_buf;
        dfe_symbols(n) = y;

        % QPSK硬判定
        d = sign(real(y)) + 1i * sign(imag(y));

        % 誤差信号（トレーニング期間は既知シンボル、以降は判定値）
        if n <= Ntrain
            e      = symbols(n) - y;
            fb_new = symbols(n);
        else
            e      = d - y;
            fb_new = d;
        end

        % LMS係数更新
        w_ff = w_ff + mu_dfe * conj(e) * ff_buf;
        w_fb = w_fb - mu_dfe * conj(e) * fb_buf;

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

end % ◀ 音源ループの終わり

%% 選択ダイバーシティ
[selected_BER,  best_rx_ber_idx]  = min(BER_matrix, [], 1);
[selected_OSNR, best_rx_osnr_idx] = max(OSNR_matrix, [], 1);

%% エリア④：最終結果の出力
clc;
disp('==========================================================================');
disp(' 比較：【選択ダイバーシティ】 vs 【PTR (16ch合成)】 vs 【PTR + DFE】');
disp('==========================================================================');
for idx_sd = 1:4
    fprintf('【音源 %d】\n', idx_sd);
    fprintf('  ➔ 選択ダイバーシティ (Best 1ch性能):\n');
    fprintf('     - 最良BER : %f (受信機 %d)\n',      selected_BER(idx_sd),  best_rx_ber_idx(idx_sd));
    fprintf('     - 最高OSNR: %.2f dB (受信機 %d)\n', selected_OSNR(idx_sd), best_rx_osnr_idx(idx_sd));
    fprintf('  ➔ PTR（16chまるごとアレイ合成):\n');
    fprintf('     - 合成BER : %f\n',    PTR_BER_vector(idx_sd));
    fprintf('     - 合成OSNR: %.2f dB\n', PTR_OSNR_vector(idx_sd));
    fprintf('  ➔ PTR + DFE (Nf=%d, Nb=%d, Ntrain=%d):\n', Nf, Nb, Ntrain);
    fprintf('     - 合成BER : %f\n',    PTRDFE_BER_vector(idx_sd));
    fprintf('     - 合成OSNR: %.2f dB\n\n', PTRDFE_OSNR_vector(idx_sd));
end
