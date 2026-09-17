%% ============================================================
%  ptr_nch_psr_sweep.m
%  「受信機数 No.Ch を増やすとPSR(集束品質)がどれだけ上がるか」を測る。
%   ・アレイ開口は保ったまま間引いて No.Ch を変える（Kida Fig.6/11の方式）
%   ・各 No.Ch で無雑音q関数のPSRを算出（矩形パルス, 現Rs）
%   ・理想の array gain 線 [1chPSR + 10log10(Nch)] と比較
%     - 実測がこの線に沿って上がる → 受信機は無相関でダイバーシティ機能
%     - 早々に飽和 → 受信機が相関（似た多重波）→ 集束が頭打ち＝真の天井
% ============================================================
arr_file = 'JamTank_Arr.arr';     % ★199ch環境のARRファイル名に変更
[Arr, Pos] = read_arrivals_asc(arr_file);

Nrr=size(Arr,1); Nrz=size(Arr,2); Nsd=size(Arr,3); NUM_RX=Nrr*Nrz;
fprintf('ARR: 距離%d × 深度%d = %dch, 音源%d\n', Nrr,Nrz,NUM_RX,Nsd);

ch_rr=zeros(NUM_RX,1); ch_rz=zeros(NUM_RX,1); cc=0;
for irr=1:Nrr, for irz=1:Nrz, cc=cc+1; ch_rr(cc)=irr; ch_rz(cc)=irz; end, end

%% ==== 設定 ====
idx_sd = 1;
Fs     = 80000;
Rs     = 5000;                    % PSRを測るシンボルレート
Lsym   = 200;
NchList = unique(min(NUM_RX, [1 2 4 6 11 20 41 60 100 150 NUM_RX]));  % 評価する受信機数

Sps=round(Fs/Rs);

%% ==== 全ch共通の最小遅延・バッファ長 ====
gmin=inf; maxd=0;
for ch=1:NUM_RX
    d=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay); d=d(d>0);
    if ~isempty(d), gmin=min(gmin,min(d)); end
end
for ch=1:NUM_RX
    d=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay); d=d(d>0)-gmin; d=d(d>=0);
    if ~isempty(d), maxd=max(maxd,max(d)); end
end
Lh=round(maxd*Fs)+8;

%% ==== 各chのh（自己相関）を事前計算 ====
qcell=cell(NUM_RX,1);
for ch=1:NUM_RX
    amp=Arr(ch_rr(ch),ch_rz(ch),idx_sd).A;
    dl =real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay);
    v=(dl>0); amp=amp(v); dl=dl(v); dl=dl-gmin;
    v=(dl>=0); amp=amp(v); dl=dl(v);
    if isempty(dl), qcell{ch}=zeros(2*Lh-1,1); continue; end
    ds=round(dl*Fs)+1; h=zeros(max(ds),1); h(ds)=amp;
    qac=conv(h,conj(flipud(h)));
    qq=zeros(2*Lh-1,1); L=min(length(qac),2*Lh-1); qq(1:L)=qac(1:L);
    qcell{ch}=qq;
end

%% ==== 単一chの平均PSR（参考: ダイバーシティ利得の基準）====
p1=zeros(NUM_RX,1);
for ch=1:NUM_RX, p1(ch)=psr_of(qcell{ch},Sps,Lsym); end
PSR_1ch = 10*log10(mean(10.^(p1/10)));   % 線形平均
fprintf('\n単一ch平均PSR = %.2f dB（Rs=%d）\n', PSR_1ch, Rs);

%% ==== No.Ch を増やしてPSR算出（開口保持で間引き）====
PSR=zeros(numel(NchList),1);
fprintf('\n%8s  %10s  %14s\n','No.Ch','PSR[dB]','理想(1ch+Ga)');
for i=1:numel(NchList)
    Nc=NchList(i);
    sel=round(linspace(1,NUM_RX,Nc));      % 開口全体に均等配置で間引き
    q=zeros(2*Lh-1,1);
    for k=1:numel(sel), q=q+qcell{sel(k)}; end
    PSR(i)=psr_of(q,Sps,Lsym);
    ideal=PSR_1ch+10*log10(Nc);
    fprintf('%8d  %10.2f  %14.2f\n', Nc, PSR(i), ideal);
end

%% ==== 描画 ====
figure('Name','PSR vs No.Ch','Position',[100 100 640 480]);
semilogx(NchList,PSR,'-o','LineWidth',1.6,'DisplayName','実測PSR'); hold on;
semilogx(NchList,PSR_1ch+10*log10(NchList),'k--','LineWidth',1.2,'DisplayName','理想: 1ch + 10log_{10}(Nch)');
grid on; xlabel('受信機数 No.Ch'); ylabel('PSR [dB] = OSNR天井');
legend('Location','best');
title(sprintf('PSR vs 受信機数（Rs=%d, 音源%d）\n理想線に沿えば無相関/飽和なら相関',Rs,idx_sd));

%% ================= ローカル関数 =================
function psr = psr_of(q, Sps, Lsym)
    [~,pk]=max(abs(q)); best=-inf;
    for ph=-floor(Sps/2):floor(Sps/2)
        lags=-Lsym:Lsym; v=zeros(numel(lags),1);
        for j=1:numel(lags)
            idx=pk+ph+lags(j)*Sps;
            if idx>=1&&idx<=length(q), v(j)=q(idx); end
        end
        main=abs(v(lags==0))^2; side=sum(abs(v(lags~=0)).^2);
        p=10*log10(main/max(side,eps));
        if p>best, best=p; end
    end
    psr=best;
end
