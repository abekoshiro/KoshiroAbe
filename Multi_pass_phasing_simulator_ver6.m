%% Bellhop データの読み込み
[Arr, Pos] = read_arrivals_asc('JamTank_Arr.arr');
BER_matrix = zeros(16,4);           % 64個の素のBERを残す箱
OSNR_matrix = zeros(16,4);          % ★64個の素のOSNRを残す箱

PTR_BER_vector     = zeros(4,1);    % PTR合成後の音源ごとのBER
PTR_OSNR_vector    = zeros(4,1);    % ★PTR合成後の音源ごとのOSNRを入れる箱
PTRDFE_BER_vector  = zeros(4,1);    % PTR+DFE後の音源ごとのBER
PTRDFE_OSNR_vector = zeros(4,1);    % ★PTR+DFE後の音源ごとのOSNRを入れる箱

plot_target_rd = 1;
plot_target_sd = 1;
idx_r = 1;

%% DFEパラメータ
Nf      = 16;    % フィードフォワードフィルタのタップ数
Nb      = 8;     % フィードバックフィルタのタップ数
mu_dfe  = 0.005; % LMSステップサイズ（収束速度と安定性のトレードオフ）
% トレーニング長（既知シンボルで係数を初期化）
Ntrain  = 200;

for idx_sd = 1:4 % ○番目の音源

    %% エリア①：【音源ループの直下】（受信機ループに入る前）
    Rs = 2000; Sps = 16; Fs = Rs * Sps; numSymbols = 2000;

    bits = randi([0 1], 2 * numSymbols, 1);
    inPhase    = 2*bits(1:2:end) - 1;
    quadrature = 2*bits(2:2:end) - 1;
    symbols    = (inPhase + 1i*quadrature); % 送信複素シンボル (s_n)
    sigPower   = mean(abs(symbols).^2);     % 送信信号の平均電力 (E[|s_n|^2])

    txSignal = zeros(numSymbols * Sps, 1);
    for i = 1:numSymbols
        txSignal((i-1)*Sps + 1 : i*Sps) = symbols(i);
    end

    % PTR用アレイ合成バケツ
    total_ptr_signal = zeros(length(txSignal) * 2, 1);

    for idx_rd = 1:16 % ○番目の受信機深度

        %% エリア②：【受信機ループの中】
        amplitudes = Arr(idx_r, idx_rd, idx_sd).A;
        delays     = Arr(idx_r, idx_rd, idx_sd).delay;

        valid_idx  = (delays > 0);
        amplitudes = amplitudes(valid_idx);
        delays     = delays(valid_idx);
        delays     = delays - min(delays);

        %% 水槽マルチパスチャネルの適用
        maxDelaySamples = round(max(delays) * Fs);
        rxSignalBase = zeros(length(txSignal) + maxDelaySamples, 1);
        for p = 1:length(delays)
            delaySamples = round(delays(p) * Fs);
            idx = delaySamples + (1:length(txSignal));
            rxSignalBase(idx) = rxSignalBase(idx) + amplitudes(p) * txSignal;
        end

        %% ノイズの付加
        SNR_dB = 20;
        sigPower_ch = mean(abs(rxSignalBase).^2);
        noisePower  = sigPower_ch / (10^(SNR_dB / 10));
        noise       = sqrt(noisePower / 2) * (randn(size(rxSignalBase)) + 1i*randn(size(rxSignalBase)));
        rxSignal    = rxSignalBase + noise;

        %% --- A. 通常の受信処理（素のBER & OSNRの計算） ---
        sampleOffset = round(Sps / 2);
        rxSymbols    = rxSignal(sampleOffset : Sps : sampleOffset + Sps*(numSymbols-1));

        alpha_single    = (symbols' * rxSymbols) / (symbols' * symbols);
        rxSymbols_comp  = rxSymbols / alpha_single;

        errPower_single               = mean(abs(symbols - rxSymbols_comp).^2);
        OSNR_matrix(idx_rd, idx_sd)   = 10 * log10(sigPower / errPower_single);

        dec_bits             = zeros(2 * numSymbols, 1);
        dec_bits(1:2:end)    = real(rxSymbols_comp) > 0;
        dec_bits(2:2:end)    = imag(rxSymbols_comp) > 0;
        BER_matrix(idx_rd, idx_sd) = sum(bits ~= dec_bits) / length(bits);

        %% --- B. PTRの個別処理と足し算 ---
        delay_samples = round(delays * Fs) + 1;
        h             = zeros(max(delay_samples), 1);
        h(delay_samples) = amplitudes;
        h_ptr_filter  = conj(flipud(h));
        z_m           = conv(rxSignal, h_ptr_filter);
        total_ptr_signal(1:length(z_m)) = total_ptr_signal(1:length(z_m)) + z_m;

    end % ◀ 受信機ループの終わり

    %% エリア③：【受信機ループを抜けた直後】（PTR合成信号の復調とBER/OSNR計算）
    [~, peak_idx] = max(abs(conv(h, conj(flipud(h)))));
    ptr_offset     = peak_idx + round(Sps / 2) - 1;

    ptrSymbols = total_ptr_signal(ptr_offset : Sps : ptr_offset + Sps*(numSymbols-1));

    alpha_ptr       = (symbols' * ptrSymbols) / (symbols' * symbols);
    ptrSymbols_comp = ptrSymbols / alpha_ptr;

    errPower_ptr              = mean(abs(symbols - ptrSymbols_comp).^2);
    PTR_OSNR_vector(idx_sd)   = 10 * log10(sigPower / errPower_ptr);

    dec_bits_ptr          = zeros(2 * numSymbols, 1);
    dec_bits_ptr(1:2:end) = real(ptrSymbols_comp) > 0;
    dec_bits_ptr(2:2:end) = imag(ptrSymbols_comp) > 0;
    PTR_BER_vector(idx_sd) = sum(bits ~= dec_bits_ptr) / length(bits);

    %% --- C. PTR + DFE処理 ---
    % PTR合成後の信号をシンボルレートでサンプリングし直す
    % （等化器への入力としてptrSymbols_compを使用：位相・振幅補正済み）
    dfe_input = ptrSymbols_comp;  % 長さ numSymbols の複素ベクトル

    % LMS-DFE フィルタ係数の初期化
    w_ff      = zeros(Nf, 1);          % フィードフォワード係数
    w_ff(1)   = 1;                     % 初期値：最初のタップを1（遅延なし応答）
    w_fb      = zeros(Nb, 1);          % フィードバック係数

    % バッファの初期化
    ff_buf    = zeros(Nf, 1);          % フィードフォワード入力バッファ
    fb_buf    = zeros(Nb, 1);          % フィードバック（判定値）バッファ

    dfe_symbols = zeros(numSymbols, 1);

    for n = 1:numSymbols
        % フィードフォワードバッファを1シンボル進める
        ff_buf = [dfe_input(n); ff_buf(1:end-1)];

        % DFE出力 = FFF出力 - FBF出力（残留ISIをキャンセル）
        y = w_ff' * ff_buf - w_fb' * fb_buf;
        dfe_symbols(n) = y;

        % QPSK硬判定
        d_real = sign(real(y));
        d_imag = sign(imag(y));
        d      = (d_real + 1i * d_imag); % 判定シンボル（±1±j）

        % LMS誤差信号の計算
        if n <= Ntrain
            % トレーニング期間：既知送信シンボルを参照
            e = symbols(n) - y;
        else
            % 判定帰還期間：硬判定結果を参照
            e = d - y;
        end

        % LMS係数更新
        w_ff = w_ff + mu_dfe * conj(e) * ff_buf;
        w_fb = w_fb - mu_dfe * conj(e) * fb_buf;

        % フィードバックバッファを更新（トレーニング中は既知シンボルで更新）
        if n <= Ntrain
            fb_buf = [symbols(n); fb_buf(1:end-1)];
        else
            fb_buf = [d; fb_buf(1:end-1)];
        end
    end

    % DFE出力のOSNR・BER計算（複素ゲイン補正）
    alpha_dfe        = (symbols' * dfe_symbols) / (symbols' * symbols);
    dfe_symbols_comp = dfe_symbols / alpha_dfe;

    errPower_dfe               = mean(abs(symbols - dfe_symbols_comp).^2);
    PTRDFE_OSNR_vector(idx_sd) = 10 * log10(sigPower / errPower_dfe);

    dec_bits_dfe          = zeros(2 * numSymbols, 1);
    dec_bits_dfe(1:2:end) = real(dfe_symbols_comp) > 0;
    dec_bits_dfe(2:2:end) = imag(dfe_symbols_comp) > 0;
    PTRDFE_BER_vector(idx_sd) = sum(bits ~= dec_bits_dfe) / length(bits);

end % ◀ 音源ループの終わり

%% 【選択ダイバーシティ】の自動適用（素の64個から最大OSNR / 最小BERを抽出）
[selected_BER,  best_rx_ber_idx]  = min(BER_matrix, [], 1);
[selected_OSNR, best_rx_osnr_idx] = max(OSNR_matrix, [], 1);

%% エリア④：最終結果の出力（BERとOSNRの比較）
clc;
disp('==========================================================================');
disp(' 比較：【選択ダイバーシティ】 vs 【PTR (16ch合成)】 vs 【PTR + DFE】');
disp('==========================================================================');
for idx_sd = 1:4
    fprintf('【音源 %d】\n', idx_sd);
    fprintf('  ➔ 選択ダイバーシティ (Best 1ch性能):\n');
    fprintf('     - 最良BER : %f (受信機 %d)\n',   selected_BER(idx_sd),  best_rx_ber_idx(idx_sd));
    fprintf('     - 最高OSNR: %.2f dB (受信機 %d)\n', selected_OSNR(idx_sd), best_rx_osnr_idx(idx_sd));
    fprintf('  ➔ PTR（16chまるごとアレイ合成):\n');
    fprintf('     - 合成BER : %f\n',    PTR_BER_vector(idx_sd));
    fprintf('     - 合成OSNR: %.2f dB\n', PTR_OSNR_vector(idx_sd));
    fprintf('  ➔ PTR + DFE (Nf=%d, Nb=%d, Ntrain=%d):\n', Nf, Nb, Ntrain);
    fprintf('     - 合成BER : %f\n',    PTRDFE_BER_vector(idx_sd));
    fprintf('     - 合成OSNR: %.2f dB\n\n', PTRDFE_OSNR_vector(idx_sd));
end
