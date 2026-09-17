%% ============================================================
%  ptr_quality_diagnostic.m  (rev2: 矩形 vs RRC 比較)
%  「PTRの集束品質」を無雑音のq関数で測り、パルス整形の効果を切り分ける。
%   ・q_delta = Σ_m conv(h_m, conj(flip(h_m)))（全受信機, 雑音なし, デルタ列）
%   ・実効シンボル応答 = q_delta ⊛ p_eff（p_eff=パルス自己相関/整形）
%       - 矩形（ver3現状）: p_eff = 幅Spsの矩形
%       - RRC整合（Kida流）: p_eff = レイズドコサイン = RRC⊛RRC（ナイキスト）
%   ・PSR(Rs) = 10log10(|q(0)|^2 / Σ_{k≠0}|q(kTs)|^2) = 無雑音OSNR上限
%   ・RRCで大きく跳ね上がれば、天井の主因は「パルス形状」だったと確定
% ============================================================
arr_file = 'JamTank_Arr.arr';     % ★199ch環境のARRファイル名に変更
[Arr, Pos] = read_arrivals_asc(arr_file);

Nrr=size(Arr,1); Nrz=size(Arr,2); Nsd=size(Arr,3); NUM_RX=Nrr*Nrz;
fprintf('ARR: 距離%d × 深度%d = %dch, 音源%d\n', Nrr,Nrz,NUM_RX,Nsd);

ch_rr=zeros(NUM_RX,1); ch_rz=zeros(NUM_RX,1); cc=0;
for irr=1:Nrr, for irz=1:Nrz, cc=cc+1; ch_rr(cc)=irr; ch_rz(cc)=irz; end, end

%% ==== 設定 ====
idx_sd  = 1;                        % 評価する音源
Fs      = 80000;                    % q関数構築の標本化周波数[Hz]
RsList  = [500 1000 2000 5000 8000];% PSRを測るシンボルレート群 [baud]
Lsym    = 200;                      % PSR計算のシンボルlag範囲(±)
beta    = 0.5;  span = 8;           % RRCパラメータ

%% ==== 無雑音q関数（デルタ列）を構築 ====
gmin=inf;
for ch=1:NUM_RX
    d=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay); d=d(d>0);
    if ~isempty(d), gmin=min(gmin,min(d)); end
end
maxd=0;
for ch=1:NUM_RX
    d=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay); d=d(d>0)-gmin; d=d(d>=0);
    if ~isempty(d), maxd=max(maxd,max(d)); end
end
Lh=round(maxd*Fs)+8; qd=zeros(2*Lh-1,1);
for ch=1:NUM_RX
    amp=Arr(ch_rr(ch),ch_rz(ch),idx_sd).A;
    dl =real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay);
    v=(dl>0); amp=amp(v); dl=dl(v); dl=dl-gmin;
    v=(dl>=0); amp=amp(v); dl=dl(v);
    if isempty(dl), continue; end
    ds=round(dl*Fs)+1; h=zeros(max(ds),1); h(ds)=amp;
    qac=conv(h,conj(flipud(h)));
    L=min(length(qac),length(qd)); qd(1:L)=qd(1:L)+qac(1:L);
end

%% ==== 矩形 vs RRC で PSR を Rs ごとに算出 ====
PSR_rect=zeros(numel(RsList),1);
PSR_rrc =zeros(numel(RsList),1);
fprintf('\n=== 無雑音 PSR（=OSNR天井）: 矩形 vs RRC ===\n');
fprintf('%8s  %10s  %10s  %8s\n','Rs','矩形[dB]','RRC[dB]','改善');
for i=1:numel(RsList)
    Sps_eff=round(Fs/RsList(i));
    % 矩形パルス（幅Sps_eff）
    p_rect=ones(Sps_eff,1);
    PSR_rect(i)=psr_shaped(qd,p_rect,Sps_eff,Lsym);
    % RRC整合（レイズドコサイン = RRC⊛RRC）
    rrcF=rrc_filter(beta,span,Sps_eff); p_rc=conv(rrcF,rrcF);
    PSR_rrc(i)=psr_shaped(qd,p_rc,Sps_eff,Lsym);
    fprintf('%8d  %10.2f  %10.2f  %+7.2f\n', RsList(i), PSR_rect(i), PSR_rrc(i), PSR_rrc(i)-PSR_rect(i));
end
fprintf('→ RRCで大きく跳ね上がれば、天井の主因は「矩形パルス」だったと確定。\n');
fprintf('  跳ね上がらなければ、集束（空間ダイバーシティ）そのものの問題。\n\n');

%% ==== 描画 ====
figure('Name','PTR集束品質: 矩形 vs RRC','Position',[80 80 1000 420]);

% (1) PSR vs Rs（矩形 vs RRC）
subplot(1,3,1);
plot(RsList,PSR_rect,'-o','LineWidth',1.5,'DisplayName','矩形(現状)'); hold on;
plot(RsList,PSR_rrc,'-s','LineWidth',1.5,'DisplayName','RRC整合(Kida流)');
grid on; xlabel('Rs [baud]'); ylabel('PSR [dB] = OSNR天井'); legend('Location','best');
title('パルス形状によるOSNR天井の違い');

% (2) 実効シンボル応答（Rs=5000, 矩形 vs RRC を重ねる）
subplot(1,3,2);
Rs0=5000; Sps0=round(Fs/Rs0);
qr=conv(qd,ones(Sps0,1)); rrcF=rrc_filter(beta,span,Sps0); qc=conv(qd,conv(rrcF,rrcF));
[qs_r,ph_r]=sample_shaped(qr,Sps0,40); [qs_c,ph_c]=sample_shaped(qc,Sps0,40);
lags=-40:40;
stem(lags-0.15,20*log10(abs(qs_r)+eps),'filled','MarkerSize',3,'DisplayName','矩形'); hold on;
stem(lags+0.15,20*log10(abs(qs_c)+eps),'filled','MarkerSize',3,'DisplayName','RRC');
grid on; ylim([-50 2]); xlabel('lag [シンボル]'); ylabel('|q(kTs)| [dB]');
legend('Location','best'); title(sprintf('Rs=%d のシンボル間隔応答',Rs0));

% (3) 参考: q_delta（時間, dB）
subplot(1,3,3);
[qpk,pk]=max(abs(qd)); t=((0:length(qd)-1)-(pk-1))/Fs*1e3;
plot(t,20*log10(abs(qd)/qpk+eps)); grid on; xlim([-5 5]); ylim([-60 2]);
xlabel('焦点からの時間 [ms]'); ylabel('|q| [dB]');
title(sprintf('q関数（%dch, 雑音なし）',NUM_RX));

sgtitle(sprintf('PTR集束品質 矩形 vs RRC（音源%d, %dch）',idx_sd,NUM_RX),'FontWeight','bold');

%% ================= ローカル関数 =================
function psr = psr_shaped(qd, pulse, Sps_eff, Lsym)
% q_delta にパルスを畳み込み、最良サンプリング位相でPSRを算出
    qs=conv(qd,pulse); [~,pk]=max(abs(qs)); best=-inf;
    for ph=-floor(Sps_eff/2):floor(Sps_eff/2)
        lags=-Lsym:Lsym; v=zeros(numel(lags),1);
        for j=1:numel(lags)
            idx=pk+ph+lags(j)*Sps_eff;
            if idx>=1&&idx<=length(qs), v(j)=qs(idx); end
        end
        main=abs(v(lags==0))^2; side=sum(abs(v(lags~=0)).^2);
        p=10*log10(main/max(side,eps));
        if p>best, best=p; end
    end
    psr=best;
end

function [v,ph_best] = sample_shaped(qs, Sps_eff, Lshow)
% 描画用：最良位相でシンボル間隔サンプルを返す（ピーク正規化）
    [qpk,pk]=max(abs(qs)); best=-inf; v=[]; ph_best=0;
    for ph=-floor(Sps_eff/2):floor(Sps_eff/2)
        idx0=pk+ph;
        if idx0-Lshow*Sps_eff<1||idx0+Lshow*Sps_eff>length(qs), continue; end
        if abs(qs(idx0))>best, best=abs(qs(idx0)); ph_best=ph; end
    end
    lags=-Lshow:Lshow; v=zeros(numel(lags),1);
    for j=1:numel(lags)
        idx=pk+ph_best+lags(j)*Sps_eff;
        if idx>=1&&idx<=length(qs), v(j)=qs(idx)/qpk; end
    end
end

function h = rrc_filter(beta, span, sps)
    N=span*sps; t=(-N/2:N/2)/sps; h=zeros(size(t));
    for i=1:numel(t)
        tt=t(i);
        if tt==0
            h(i)=1-beta+4*beta/pi;
        elseif abs(abs(tt)-1/(4*beta))<1e-8
            h(i)=(beta/sqrt(2))*((1+2/pi)*sin(pi/(4*beta))+(1-2/pi)*cos(pi/(4*beta)));
        else
            num=sin(pi*tt*(1-beta))+4*beta*tt.*cos(pi*tt*(1+beta));
            den=pi*tt.*(1-(4*beta*tt).^2); h(i)=num/den;
        end
    end
    h=h/sqrt(sum(h.^2));
end
