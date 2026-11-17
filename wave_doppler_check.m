%% ============================================================
%  wave_doppler_check.m
%  時変チャネル化で「8枚のスナップショットが足りているか」を検証する診断。
%
%  問題意識:
%   水面が動くと海面反射パスの経路長が音響波長オーダーで変化し、
%   複素CIRの位相が速く回転する。隣接スナップショット間の位相回転が
%   π を超えると、複素CIRの線形補間はエイリアシング（偽の振幅ノッチ）に
%   なり、時変チャネルが正しく再現できない。
%
%  やること:
%   各受信機で隣接スナップショットのCIR相関 c_k=<h_k,h_{k+1}>/(|h_k||h_{k+1}|)
%   を計算し、
%     ・位相ステップ |angle(c_k)| … これが π に近い/超えるとアウト
%     ・相関振幅 |c_k|          … 1に近いほど滑らか、低いと急変
%   を全受信機で集計。推奨スナップショット数も提示する。
%
%  依存: read_arrivals_asc.m
% ============================================================
arr_base = 'JamTank_Arr';
n_snap   = 8;
idx_sd   = 1;          % 評価する音源
Fs       = 80000;      % Wave_TimeVarying_simulator と同じサンプリング
Lfrac    = 16;

%% 読込
Arr_all=cell(n_snap,1);
for k=1:n_snap
    [Arr_all{k},~]=read_arrivals_asc(sprintf('%s_wave%02d.arr',arr_base,k));
end
Nrr=size(Arr_all{1},1); Nrz=size(Arr_all{1},2); NUM_RX=Nrr*Nrz;
ch_rr=zeros(NUM_RX,1); ch_rz=zeros(NUM_RX,1); cc=0;
for irr=1:Nrr, for irz=1:Nrz, cc=cc+1; ch_rr(cc)=irr; ch_rz(cc)=irz; end, end

%% 共通CIR長
gmin=inf; Lmax=0;
for k=1:n_snap, for ch=1:NUM_RX
    d=real(Arr_all{k}(ch_rr(ch),ch_rz(ch),idx_sd).delay); d=d(d>0);
    if ~isempty(d), gmin=min(gmin,min(d)); end
end, end
for k=1:n_snap, for ch=1:NUM_RX
    d=real(Arr_all{k}(ch_rr(ch),ch_rz(ch),idx_sd).delay); d=d(d>0)-gmin; d=d(d>=0);
    if ~isempty(d), Lmax=max(Lmax,round(max(d)*Fs)+1); end
end, end
Lh=Lmax+2*Lfrac+4;

%% 各受信機: スナップショット列のCIRを作り、隣接相関を計算
step_ph=[]; corr_mag=[];   % 全(受信機×ステップ)のサンプルを貯める
per_rx_maxstep=nan(NUM_RX,1);
for ch=1:NUM_RX
    H=zeros(Lh,n_snap);
    for k=1:n_snap
        amp=Arr_all{k}(ch_rr(ch),ch_rz(ch),idx_sd).A;
        dl =real(Arr_all{k}(ch_rr(ch),ch_rz(ch),idx_sd).delay);
        v=(dl>0); amp=amp(v); dl=dl(v); dl=dl-gmin; v=(dl>=0); amp=amp(v); dl=dl(v);
        if isempty(dl), continue; end
        H(:,k)=build_cir_frac(amp,dl*Fs+Lfrac,Lh,Lfrac);
    end
    if all(H(:)==0), continue; end
    % 周期的に隣接相関（k=n_snap と k=1 も隣接）
    mx=0;
    for k=1:n_snap
        a=H(:,k); b=H(:,mod(k,n_snap)+1);
        na=norm(a); nb=norm(b); if na<eps||nb<eps, continue; end
        c=(a'*b)/(na*nb);
        step_ph(end+1)=abs(angle(c)); %#ok<SAGROW>
        corr_mag(end+1)=abs(c);        %#ok<SAGROW>
        mx=max(mx,abs(angle(c)));
    end
    per_rx_maxstep(ch)=mx;
end

%% 集計・判定
step_ph=step_ph(:); corr_mag=corr_mag(:);
med_ph=median(step_ph); max_ph=max(step_ph);
med_cm=median(corr_mag);
frac_over=mean(step_ph>pi);          % 位相ステップが π を超えた割合
% 推奨: 位相ステップを pi/4 以下に抑えるのに必要なスナップショット数
rec_nsnap=ceil(n_snap*max_ph/(pi/4));

fprintf('==================================================================\n');
fprintf(' 波スナップショット十分性チェック（音源%d, %dch, %d枚/周期）\n',idx_sd,NUM_RX,n_snap);
fprintf('==================================================================\n');
fprintf(' 隣接スナップ間 位相ステップ |Δφ|:  中央値 %.2f rad (%.0f°) / 最大 %.2f rad (%.0f°)\n',...
        med_ph,rad2deg(med_ph),max_ph,rad2deg(max_ph));
fprintf(' 隣接スナップ間 CIR相関 |c|     :  中央値 %.3f（1=滑らか, 低=急変）\n',med_cm);
fprintf(' 位相ステップが π を超えた割合  :  %.0f%%\n',frac_over*100);
fprintf('------------------------------------------------------------------\n');
if max_ph < pi/2
    fprintf(' 判定: ○ 8枚で十分に滑らか（補間は信頼できる）。\n');
    fprintf('       → 時変の劣化は【物理的に本物】。凍結DFEが波に弱いことの表れ。\n');
elseif max_ph < pi
    fprintf(' 判定: △ やや粗い。最大 %.0f° 回転。補間に誤差が混じる可能性。\n',rad2deg(max_ph));
    fprintf('       → n_snap=%d 程度に増やして再確認を推奨。\n',rec_nsnap);
else
    fprintf(' 判定: × 不足（最大 %.0f°>180°）。複素CIR補間がエイリアシング。\n',rad2deg(max_ph));
    fprintf('       → 現在の 4dB崩壊はサンプリング不足のアーティファクトの疑い大。\n');
    fprintf('       → make_wave_ati の n_snap を %d 以上にして Bellhop を回し直す。\n',rec_nsnap);
end
fprintf('==================================================================\n');

%% 図
figure('Name','波スナップショット十分性','Position',[100 100 900 360]);
subplot(1,2,1);
histogram(rad2deg(step_ph),20); hold on;
yl=ylim; plot([180 180],yl,'r--','LineWidth',1.5);
text(182,yl(2)*0.9,'π (アウト境界)','Color','r');
xlabel('隣接スナップ間 位相ステップ |Δφ| [deg]'); ylabel('件数');
title(sprintf('位相ステップ分布（中央%.0f° 最大%.0f°）',rad2deg(med_ph),rad2deg(max_ph)));
grid on;
subplot(1,2,2);
histogram(corr_mag,20);
xlabel('隣接スナップ間 CIR相関 |c|'); ylabel('件数');
title(sprintf('CIR相関分布（中央%.3f）',med_cm)); grid on;
sgtitle('波スナップショットが時変チャネルを解像できているかの診断','FontWeight','bold');

%% ローカル関数
function h = build_cir_frac(amp, cs, Lh, L)
    h=zeros(Lh,1);
    for p=1:numel(amp)
        c=cs(p); m=(floor(c)-L:floor(c)+L).';
        ok=m>=0 & m<=Lh-1; m=m(ok); if isempty(m), continue; end
        x=m-c; w=0.5*(1+cos(pi*x/(L+1)));
        h(m+1)=h(m+1)+amp(p)*(sinc(x).*w);
    end
end
