%% ============================================================
%  PTR 診断スクリプト（音源1・1音源のみ／DFEトグル付き）
%  目的：q関数の鋭さ・PTR出力コンスタレーションを「目で見る」
% ============================================================
[Arr, Pos] = read_arrivals_asc('JamTank_Arr.arr');

idx_r  = 1;
idx_sd = 1;          % 診断する音源番号
USE_DFE = false;     % ★まずはDFE無し(false)でPTR単体を確認する

%% 送信信号
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
for k = 1:16
    d = Arr(idx_r,k,idx_sd).delay; d = d(d>0);
    if ~isempty(d), global_min_delay = min(global_min_delay, min(d)); end
end

ptr_buf_len      = length(txSignal) * 4;
total_ptr_signal = zeros(ptr_buf_len,1);
total_autocorr   = zeros(ptr_buf_len,1);

for idx_rd = 1:16
    amplitudes = Arr(idx_r,idx_rd,idx_sd).A;
    delays     = Arr(idx_r,idx_rd,idx_sd).delay;
    v = (delays>0); amplitudes=amplitudes(v); delays=delays(v);
    delays = delays - global_min_delay;
    v = (delays>=0); amplitudes=amplitudes(v); delays=delays(v);
    if isempty(delays), continue; end

    maxDelaySamples = round(max(delays)*Fs);
    rxBase = zeros(length(txSignal)+maxDelaySamples,1);
    for p=1:length(delays)
        ds = round(delays(p)*Fs);
        rxBase(ds+(1:length(txSignal))) = rxBase(ds+(1:length(txSignal))) + amplitudes(p)*txSignal;
    end
    SNR_dB=20; sp=mean(abs(rxBase).^2); np=sp/10^(SNR_dB/10);
    rx = rxBase + sqrt(np/2)*(randn(size(rxBase))+1i*randn(size(rxBase)));

    ds2 = round(delays*Fs)+1;
    h = zeros(max(ds2),1); h(ds2)=amplitudes;
    g = conj(flipud(h));
    z = conv(rx,g);
    L = min(length(z),ptr_buf_len);
    total_ptr_signal(1:L) = total_ptr_signal(1:L) + z(1:L);

    hac = conv(h,g);
    La  = min(length(hac),ptr_buf_len);
    total_autocorr(1:La) = total_autocorr(1:La) + hac(1:La);
end

%% q関数（合成自己相関）の評価
[qpeak, peak_idx] = max(abs(total_autocorr));
% シンボル間隔(Sps)でのサイドローブを抽出
lags      = (-numSymbols:numSymbols);
q_at_sym  = zeros(length(lags),1);
for ii=1:length(lags)
    idx = peak_idx + lags(ii)*Sps;
    if idx>=1 && idx<=length(total_autocorr)
        q_at_sym(ii) = total_autocorr(idx);
    end
end
mainlobe   = abs(q_at_sym(lags==0))^2;
sidelobe   = sum(abs(q_at_sym(lags~=0)).^2);
PSR_dB     = 10*log10(mainlobe / sidelobe);   % ピーク対サイドローブ比（理論OSNR上限の目安）

fprintf('q関数 ピーク対サイドローブ比 (symbol-spaced): %.2f dB\n', PSR_dB);
fprintf(' → これがOSNRのおおよその上限です。低ければPTRの焦点が甘い証拠。\n\n');

%% PTRシンボル抽出
ptr_offset = peak_idx + round(Sps/2) - 1;
ptrSymbols = total_ptr_signal(ptr_offset : Sps : ptr_offset + Sps*(numSymbols-1));
alpha_ptr  = (symbols'*ptrSymbols)/(symbols'*symbols);
ptrComp    = ptrSymbols/alpha_ptr;
osnr_ptr   = 10*log10(sigPower/mean(abs(symbols-ptrComp).^2));
fprintf('PTR単体 OSNR = %.2f dB\n', osnr_ptr);

%% 図1：q関数（シンボル間隔サンプリング）
figure('Name','q-function');
stem(lags, abs(q_at_sym)/qpeak, 'filled'); xlim([-100 100]);
xlabel('シンボル遅延'); ylabel('|q| / |q_{peak}|');
title(sprintf('合成q関数（PSR=%.1f dB）：理想は中央のみ1で他が0', PSR_dB)); grid on;

%% 図2：PTR出力コンスタレーション
figure('Name','PTR constellation');
plot(real(ptrComp), imag(ptrComp), '.'); axis equal; grid on;
xlabel('I'); ylabel('Q'); title(sprintf('PTR出力（OSNR=%.1f dB）：4点に集まれば成功', osnr_ptr));

%% （任意）DFE
if USE_DFE
    Nf=80; Nb=80; mu=0.001; Ntr=500;
    wff=zeros(Nf,1); wff(1)=1; wfb=zeros(Nb,1);
    fbuf=zeros(Nf,1); bbuf=zeros(Nb,1); out=zeros(numSymbols,1);
    for n=1:numSymbols
        fbuf=[ptrSymbols(n);fbuf(1:end-1)];
        y=wff'*fbuf - wfb'*bbuf; out(n)=y;
        d=sign(real(y))+1i*sign(imag(y));
        if n<=Ntr, e=symbols(n)-y; nb=symbols(n); else, e=d-y; nb=d; end
        wff=wff+mu*conj(e)*fbuf; wfb=wfb-mu*conj(e)*bbuf;
        bbuf=[nb;bbuf(1:end-1)];
    end
    a=(symbols'*out)/(symbols'*symbols); oc=out/a;
    osnr_dfe=10*log10(sigPower/mean(abs(symbols-oc).^2));
    fprintf('PTR+DFE OSNR = %.2f dB\n', osnr_dfe);
    figure('Name','PTR+DFE constellation');
    plot(real(oc),imag(oc),'.'); axis equal; grid on;
    title(sprintf('PTR+DFE出力（OSNR=%.1f dB）', osnr_dfe));
end
