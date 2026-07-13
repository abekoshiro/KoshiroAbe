%% ============================================================
%  input_signal_diagnostic.m
%  「入力側」の検算：変調前のチャネル情報と送受信のやりとりを可視化
%  出力内容
%    ① CIR（チャネルインパルス応答）        … 想定パスが入っているか
%    ② 遅延広がり（最大超過遅延 / RMS）      … ISIの大きさ
%    ③ コヒーレンス帯域幅 Bc                 … フラット/周波数選択性の判定
%    ④ 周波数応答 H(f)（CIRのFFT）           … フェージングのノッチ
%    ⑤ パワー遅延プロファイル PDP
%    ⑥ 送信信号のスペクトル・占有帯域
%    ⑦ 1対1リンクの畳込み可視化（tx→ch→rx）
% ============================================================
[Arr, Pos] = read_arrivals_asc('JamTank_Arr.arr');

Nrr = size(Arr,1); Nrz = size(Arr,2); Nsd = size(Arr,3);
NUM_RX = Nrr*Nrz;

%% ==== 確認対象を選択（ここを変えて他の音源/受信機も確認）====
idx_sd = 1;   % 音源番号
idx_rr = 1;   % 受信距離インデックス（1..Nrr）
idx_rz = 1;   % 受信深度インデックス（1..Nrz）

%% ==== 送信信号パラメータ（simulatorと同じ）====
Rs = 500; Sps = 16; Fs = Rs*Sps; numSymbols = 2000;
Ts = 1/Rs;                    % シンボル長 [s]

fprintf('=== 入力診断: 音源%d, 受信(距離idx%d, 深度idx%d) ===\n', idx_sd, idx_rr, idx_rz);
if isfield(Pos,'r')&&isfield(Pos.r,'r'), fprintf('  受信距離: %.1f m\n', Pos.r.r(idx_rr)); end
if isfield(Pos,'r')&&isfield(Pos.r,'z'), fprintf('  受信深度: %.2f m\n', Pos.r.z(idx_rz)); end
fprintf('  Rs=%d baud, Sps=%d, Fs=%d Hz, シンボル長Ts=%.3f ms\n\n', Rs, Sps, Fs, Ts*1e3);

%% ==== ① CIR の構築（変調前・Arrそのまま）====
amp = Arr(idx_rr, idx_rz, idx_sd).A;
dl  = real(Arr(idx_rr, idx_rz, idx_sd).delay);
v = (dl>0); amp=amp(v); dl=dl(v);
dl = dl - min(dl);            % 最初のパスを0基準に
[dl, si] = sort(dl); amp = amp(si);
Npath = numel(dl);

fprintf('① 到来パス数: %d 本\n', Npath);
fprintf('   No.  遅延[ms]   |振幅|      位相[deg]\n');
for p = 1:Npath
    fprintf('   %2d   %7.3f   %8.2e   %7.1f\n', p, dl(p)*1e3, abs(amp(p)), angle(amp(p))*180/pi);
end

% サンプル格子上のCIR
ds = round(dl*Fs)+1;
h  = zeros(max(ds),1);
for p=1:Npath, h(ds(p)) = h(ds(p)) + amp(p); end
t_h = (0:length(h)-1)/Fs;     % 時間軸 [s]

%% ==== ② 遅延広がり ====
tau_max = max(dl) - min(dl);                       % 最大超過遅延
Pk      = abs(amp).^2;                             % 各パスのパワー
tau_bar = sum(dl.*Pk)/sum(Pk);                     % 平均遅延（パワー重心）
tau_rms = sqrt( sum((dl-tau_bar).^2 .* Pk)/sum(Pk) ); % RMS遅延広がり
fprintf('\n② 最大超過遅延 τ_max = %.3f ms (= %.1f シンボル長)\n', tau_max*1e3, tau_max/Ts);
fprintf('   RMS遅延広がり στ  = %.3f ms (= %.1f シンボル長)\n', tau_rms*1e3, tau_rms/Ts);

%% ==== ③ コヒーレンス帯域幅 ====
Bc = 1/(2*pi*tau_rms);
sig_bw = Rs;                  % QPSK矩形の目安帯域（メインローブ片側≈Rs）
fprintf('\n③ コヒーレンス帯域幅 Bc ≈ %.1f Hz\n', Bc);
fprintf('   信号帯域(目安) ≈ %.1f Hz\n', sig_bw);
if sig_bw > Bc
    fprintf('   → 信号帯域 > Bc : 周波数選択性フェージング（ISIあり, PTR/等化が有効）\n');
else
    fprintf('   → 信号帯域 < Bc : フラットフェージング（ISI小）\n');
end

%% ==== ④ 周波数応答 H(f) ====
Nfft = 4096;
H = fftshift(fft(h, Nfft));
f  = (-Nfft/2:Nfft/2-1)*(Fs/Nfft);
Hdb = 20*log10(abs(H)/max(abs(H))+eps);

%% ==== ⑥ 送信信号（QPSK矩形）とスペクトル ====
bits = randi([0 1], 2*numSymbols,1);
symbols = (2*bits(1:2:end)-1)+1i*(2*bits(2:2:end)-1);
txSignal = zeros(numSymbols*Sps,1);
for i=1:numSymbols, txSignal((i-1)*Sps+1:i*Sps)=symbols(i); end
Ntx = 4096;
TX = fftshift(fft(txSignal(1:min(Ntx,end)), Ntx));
ftx = (-Ntx/2:Ntx/2-1)*(Fs/Ntx);
TXdb = 20*log10(abs(TX)/max(abs(TX))+eps);

%% ==== ⑦ 1対1リンク畳込み（tx→ch→rx）====
rx_clean = conv(txSignal, h);                     % チャネル通過（雑音なし）
sp = mean(abs(rx_clean).^2); np = sp/10^(20/10);  % SNR=20dB
rx_noisy = rx_clean + sqrt(np/2)*(randn(size(rx_clean))+1i*randn(size(rx_clean)));
t_tx = (0:length(txSignal)-1)/Fs;
t_rx = (0:length(rx_clean)-1)/Fs;

%% ============ 描画 ============
figure('Name','入力診断（CIR/遅延/周波数応答）','Position',[60 60 1200 760]);

% ① CIR
subplot(2,3,1);
stem(t_h*1e3, abs(h), 'filled','MarkerSize',3); grid on;
xlabel('遅延 [ms]'); ylabel('|h|'); xlim([0 max(tau_max*1e3*1.2, 0.5)]);
title(sprintf('① CIR（%d パス）', Npath));

% ⑤ PDP（dB）
subplot(2,3,2);
stem(dl*1e3, 10*log10(Pk/max(Pk)), 'filled','MarkerSize',3); grid on;
xlabel('遅延 [ms]'); ylabel('相対パワー [dB]');
title('⑤ パワー遅延プロファイル');

% ④ 周波数応答
subplot(2,3,3);
plot(f, Hdb); grid on; xlabel('周波数 [Hz]'); ylabel('|H(f)| [dB]');
xlim([-Fs/2 Fs/2]);
title(sprintf('④ 周波数応答  Bc≈%.0f Hz', Bc));
hold on;
yl=ylim;
plot([-sig_bw -sig_bw],yl,'r--'); plot([sig_bw sig_bw],yl,'r--');
legend({'|H(f)|','信号帯域±Rs'},'Location','south'); hold off;

% ⑥ 送信スペクトル
subplot(2,3,4);
plot(ftx, TXdb); grid on; xlabel('周波数 [Hz]'); ylabel('|TX(f)| [dB]');
xlim([-Fs/2 Fs/2]); title('⑥ 送信信号スペクトル（QPSK矩形）');

% ⑦ 送信波形（拡大）
subplot(2,3,5);
Nshow = 20*Sps;   % 先頭20シンボル
plot(t_tx(1:Nshow)*1e3, real(txSignal(1:Nshow)),'b'); hold on;
plot(t_tx(1:Nshow)*1e3, imag(txSignal(1:Nshow)),'r'); grid on;
xlabel('時間 [ms]'); ylabel('振幅'); legend({'I','Q'});
title('⑦a 送信波形（先頭20シンボル）');

% ⑦ 受信波形（同区間, チャネル通過後）
subplot(2,3,6);
Nshow2 = Nshow + length(h);
plot(t_rx(1:Nshow2)*1e3, real(rx_noisy(1:Nshow2)),'b'); hold on;
plot(t_rx(1:Nshow2)*1e3, real(rx_clean(1:Nshow2)),'k','LineWidth',1); grid on;
xlabel('時間 [ms]'); ylabel('Re'); legend({'雑音あり','雑音なし(tx⊛h)'});
title('⑦b 受信波形 = tx ⊛ h + 雑音');

sgtitle(sprintf('入力診断: 音源%d → 受信(距離idx%d,深度idx%d)  τ_{max}=%.2fms, στ=%.3fms', ...
        idx_sd, idx_rr, idx_rz, tau_max*1e3, tau_rms*1e3), 'FontWeight','bold');
