%% ============================================================
%  ptr_quality_diagnostic.m
%  「PTRの集束品質」を無雑音のq関数で測る。
%   ・q関数 = Σ_m conv(h_m, conj(flip(h_m)))（全受信機, 雑音なし）
%   ・PSR(Rs) = 10log10( |q(0)|^2 / Σ_{k≠0}|q(kTs)|^2 )
%       = そのシンボルレートでの focusing-limited OSNR 上限（雑音なし理論値）
%   ・Rsを掃引 → 「Rsを上げると天井が下がる」トレードオフを可視化
%   ・全ch単純PTR と 距離群別(参考)の両方を評価
% ============================================================
arr_file = 'JamTank_Arr.arr';     % ★199ch環境のARRファイル名に変更
[Arr, Pos] = read_arrivals_asc(arr_file);

Nrr=size(Arr,1); Nrz=size(Arr,2); Nsd=size(Arr,3); NUM_RX=Nrr*Nrz;
fprintf('ARR: 距離%d × 深度%d = %dch, 音源%d\n', Nrr,Nrz,NUM_RX,Nsd);

ch_rr=zeros(NUM_RX,1); ch_rz=zeros(NUM_RX,1); cc=0;
for irr=1:Nrr, for irz=1:Nrz, cc=cc+1; ch_rr(cc)=irr; ch_rz(cc)=irz; end, end

%% ==== 設定 ====
idx_sd  = 1;                    % 評価する音源
Fs      = 80000;                % q関数構築の標本化周波数[Hz]（高いほど遅延分解能良）
RsList  = [500 1000 2000 5000 8000];  % PSRを測るシンボルレート群 [baud]
Lsym    = 200;                  % PSR計算で見るシンボルlag範囲（±Lsym）

%% ==== 全受信機共通の最小遅延で整列してq関数を構築（雑音なし）====
gmin=inf;
for ch=1:NUM_RX
    d=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay); d=d(d>0);
    if ~isempty(d), gmin=min(gmin,min(d)); end
end

% 最大遅延からバッファ長を決定
maxd=0;
for ch=1:NUM_RX
    d=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay); d=d(d>0)-gmin; d=d(d>=0);
    if ~isempty(d), maxd=max(maxd,max(d)); end
end
Lh=round(maxd*Fs)+8;
qsum=zeros(2*Lh-1,1);

for ch=1:NUM_RX
    amp=Arr(ch_rr(ch),ch_rz(ch),idx_sd).A;
    dl =real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay);
    v=(dl>0); amp=amp(v); dl=dl(v); dl=dl-gmin;
    v=(dl>=0); amp=amp(v); dl=dl(v);
    if isempty(dl), continue; end
    ds=round(dl*Fs)+1; h=zeros(max(ds),1); h(ds)=amp;
    qac=conv(h,conj(flipud(h)));
    L=min(length(qac),length(qsum)); qsum(1:L)=qsum(1:L)+qac(1:L);
end

[qpk,pk]=max(abs(qsum));
qn=qsum/qpk;                     % ピーク正規化

%% ==== PSR を Rs ごとに算出 ====
PSR=zeros(numel(RsList),1);
fprintf('\n=== 無雑音 q関数 PSR（= focusing-limited OSNR上限）===\n');
for i=1:numel(RsList)
    Sps_eff=round(Fs/RsList(i));      % そのRsでのシンボル間隔[サンプル]
    lags=-Lsym:Lsym; qs=zeros(numel(lags),1);
    for j=1:numel(lags)
        idx=pk+lags(j)*Sps_eff;
        if idx>=1&&idx<=length(qsum), qs(j)=qsum(idx); end
    end
    main=abs(qs(lags==0))^2; side=sum(abs(qs(lags~=0)).^2);
    PSR(i)=10*log10(main/max(side,eps));
    fprintf('  Rs=%5d baud (シンボル間隔%3dサンプル): PSR = %6.2f dB\n', RsList(i), Sps_eff, PSR(i));
end
fprintf('→ このPSRが、その各Rsでの「雑音ゼロでも超えられないOSNR天井」です。\n');
fprintf('  Rsを上げるほどPSRが下がる＝データレートとOSNRのトレードオフ。\n\n');

%% ==== 描画 ====
figure('Name','PTR集束品質 診断','Position',[80 80 1000 420]);

% (1) q関数（時間, dB）: メインローブとpedestalを見る
subplot(1,3,1);
t=((0:length(qsum)-1)-(pk-1))/Fs*1e3;   % 焦点を0[ms]に
plot(t,20*log10(abs(qn)+eps)); grid on; xlim([-5 5]);
xlabel('焦点からの時間 [ms]'); ylabel('|q| [dB]'); ylim([-60 2]);
title(sprintf('q関数（%dch, 雑音なし）\n中央が鋭くpedestalが低いほど良',NUM_RX));

% (2) PSR vs Rs
subplot(1,3,2);
plot(RsList,PSR,'-o','LineWidth',1.5); grid on;
xlabel('シンボルレート Rs [baud]'); ylabel('PSR [dB] = OSNR天井');
title('Rsを上げると天井が下がる');

% (3) 現Rsのシンボル間隔q（stem）
subplot(1,3,3);
Rs0=5000; Sps0=round(Fs/Rs0); lags=-40:40; qs=zeros(numel(lags),1);
for j=1:numel(lags), idx=pk+lags(j)*Sps0; if idx>=1&&idx<=length(qsum), qs(j)=qsum(idx); end, end
stem(lags,20*log10(abs(qs)/qpk+eps),'filled','MarkerSize',3); grid on;
xlabel('lag [シンボル]'); ylabel('|q(kTs)| [dB]'); ylim([-50 2]);
title(sprintf('Rs=%d のシンボル間隔q\n(残留ISI構造)',Rs0));
sgtitle(sprintf('PTR集束品質（音源%d, %dch）',idx_sd,NUM_RX),'FontWeight','bold');
