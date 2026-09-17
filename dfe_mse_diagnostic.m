%% ============================================================
%  dfe_mse_diagnostic.m
%  DFEの収束をMSE学習曲線で確認する（星座図で点が飛ぶ原因の切り分け）
%    ・|e[n]|^2 をシンボル番号 n に対してプロット（学習曲線）
%    ・理想: トレーニング中にMSEが下がり、以降フラットに低い
%    ・スパイク/バースト → 誤り伝播（判定ミスがFB経由で伝播）＝点が飛ぶ原因
%    ・星座図を時間(n)で色分け → いつ飛んだかをMSEスパイクと対応づけ
% ============================================================
arr_file = 'JamTank_Arr_ver3.arr';    % ★使用中のARRに合わせる
[Arr, Pos] = read_arrivals_asc(arr_file);
Nrr=size(Arr,1); Nrz=size(Arr,2); Nsd=size(Arr,3); NUM_RX=Nrr*Nrz;

%% ==== 設定（ver3に合わせる）====
idx_sd = 1;
Rs=5000; Sps=16; Fs=Rs*Sps; numSymbols=2000;
Nf=20; Nb=15; lambda=0.999; delta=0.01; Ntrain=600;
SNR_dB=20;
win=31;                 % 学習曲線の移動平均窓[シンボル]

ch_rr=zeros(NUM_RX,1); ch_rz=zeros(NUM_RX,1); cc=0;
for irr=1:Nrr, for irz=1:Nrz, cc=cc+1; ch_rr(cc)=irr; ch_rz(cc)=irz; end, end
if Nrr>=2, NUM_SUB=Nrr; ch_sub=ch_rr;
else, NUM_SUB=2; SUB=ceil(NUM_RX/NUM_SUB); ch_sub=min(NUM_SUB,ceil((1:NUM_RX)'/SUB)); end

%% ==== 送信・チャネル・SCM-PTR（1試行）====
bits=randi([0 1],2*numSymbols,1);
symbols=(2*bits(1:2:end)-1)+1i*(2*bits(2:2:end)-1);
txSignal=zeros(numSymbols*Sps,1);
for i=1:numSymbols, txSignal((i-1)*Sps+1:i*Sps)=symbols(i); end

gmin=inf;
for ch=1:NUM_RX
    d=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay); d=d(d>0);
    if ~isempty(d), gmin=min(gmin,min(d)); end
end
ptr_buf_len=length(txSignal)*4;
sub_ptr=zeros(ptr_buf_len,NUM_SUB); sub_ac=zeros(ptr_buf_len,NUM_SUB);
for ch=1:NUM_RX
    amp=Arr(ch_rr(ch),ch_rz(ch),idx_sd).A; dl=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay);
    v=(dl>0); amp=amp(v); dl=dl(v); dl=dl-gmin; v=(dl>=0); amp=amp(v); dl=dl(v);
    if isempty(dl), continue; end
    mds=round(max(dl)*Fs); rxBase=zeros(length(txSignal)+mds,1);
    for p=1:length(dl), ds=round(dl(p)*Fs); rxBase(ds+(1:length(txSignal)))=rxBase(ds+(1:length(txSignal)))+amp(p)*txSignal; end
    sp=mean(abs(rxBase).^2); np=sp/10^(SNR_dB/10);
    rx=rxBase+sqrt(np/2)*(randn(size(rxBase))+1i*randn(size(rxBase)));
    g=conj(flipud((function_h(dl,amp,Fs)))); z=conv(rx,g); hac=conv(function_h(dl,amp,Fs),g);
    Lz=min(length(z),ptr_buf_len); La=min(length(hac),ptr_buf_len);
    si=ch_sub(ch);
    sub_ptr(1:Lz,si)=sub_ptr(1:Lz,si)+z(1:Lz); sub_ac(1:La,si)=sub_ac(1:La,si)+hac(1:La);
end
scm=zeros(numSymbols,1); wsum=0;
for s=1:NUM_SUB
    [~,pks]=max(abs(sub_ac(:,s))); offs=pks+round(Sps/2)-1;
    if offs+Sps*(numSymbols-1)>ptr_buf_len, continue; end
    ss=sub_ptr(offs:Sps:offs+Sps*(numSymbols-1),s);
    as=(symbols'*ss)/(symbols'*symbols); sc=ss/as; es=mean(abs(symbols-sc).^2); ws=1/es;
    scm=scm+ws*sc; wsum=wsum+ws;
end
scm=scm/wsum;

%% ==== RLS-DFE（誤差e[n]も記録）====
[dfe_out, err] = rls_dfe_log(scm, symbols, Nf, Nb, lambda, delta, Ntrain);
mse_inst = abs(err).^2;                          % 瞬時二乗誤差
mse_ma   = movmean(mse_inst, win);               % 移動平均MSE
ss_idx   = (Ntrain+1):numSymbols;
mse_ss   = mean(mse_inst(ss_idx));               % 定常MSE（トレーニング後）

fprintf('=== DFE MSE診断（音源%d, Rs=%d, Nf=%d Nb=%d Ntrain=%d）===\n', idx_sd,Rs,Nf,Nb,Ntrain);
fprintf('  定常MSE（トレーニング後平均）= %.4e  = %.2f dB\n', mse_ss, 10*log10(mse_ss));
% バースト（誤り伝播）検出：定常MSEの10倍を超えるシンボル
burst = find(mse_inst(ss_idx) > 10*mse_ss);
fprintf('  バースト（定常MSEの10倍超）: %d 個 / %d シンボル\n', numel(burst), numel(ss_idx));

%% ==== 描画 ====
figure('Name','DFE MSE学習曲線','Position',[80 80 1050 440]);

% (1) 学習曲線（全体, dB）
subplot(1,3,1);
plot(1:numSymbols, 10*log10(mse_inst+eps), 'Color',[0.7 0.8 1]); hold on;
plot(1:numSymbols, 10*log10(mse_ma+eps), 'b','LineWidth',1.5);
xline(Ntrain,'r--','トレーニング終了');
yline(10*log10(mse_ss),'k:','定常MSE');
grid on; xlabel('シンボル番号 n'); ylabel('MSE = |e[n]|^2 [dB]');
title('DFE学習曲線（下がって平らなら収束）');
legend({'瞬時','移動平均'},'Location','best');

% (2) 定常部の拡大（バーストを見る）
subplot(1,3,2);
plot(ss_idx, 10*log10(mse_inst(ss_idx)+eps), 'Color',[1 0.7 0.6]); hold on;
plot(ss_idx, 10*log10(mse_ma(ss_idx)+eps), 'r','LineWidth',1.3);
yline(10*log10(10*mse_ss),'k--','バースト閾値');
grid on; xlabel('シンボル番号 n'); ylabel('MSE [dB]');
title(sprintf('定常部（バースト%d個）\n頻繁なら誤り伝播で点が飛ぶ',numel(burst)));

% (3) 星座図を時間で色分け
subplot(1,3,3);
scatter(real(dfe_out(ss_idx)),imag(dfe_out(ss_idx)),18,ss_idx,'filled'); hold on;
plot([1 -1 1 -1],[1 1 -1 -1],'r+','MarkerSize',12,'LineWidth',2);
axis equal; grid on; colormap(jet); cb=colorbar; ylabel(cb,'シンボル番号 n');
lim=max(3,ceil(max(abs([real(dfe_out(ss_idx));imag(dfe_out(ss_idx))]))));
xlim([-lim lim]); ylim([-lim lim]); xlabel('I'); ylabel('Q');
title('DFE出力星座図（時間で色分け）');
sgtitle('DFE収束のMSE診断','FontWeight','bold');

%% ================= ローカル関数 =================
function h = function_h(dl, amp, Fs)
    ds=round(dl*Fs)+1; h=zeros(max(ds),1); h(ds)=amp;
end

function [out, err] = rls_dfe_log(x, symbols, Nf, Nb, lambda, delta, Ntrain)
% RLS-DFE。out=等化出力, err=各シンボルの誤差 e[n]=s_ref-y
    N=length(x); Ntap=Nf+Nb;
    w=zeros(Ntap,1); w(1)=1; P=(1/delta)*eye(Ntap);
    ff=zeros(Nf,1); fb=zeros(Nb,1); out=zeros(N,1); err=zeros(N,1);
    for n=1:N
        ff=[x(n); ff(1:end-1)]; u=[ff; -fb];
        y=w'*u; out(n)=y;
        d=sign(real(y))+1i*sign(imag(y));
        if n<=Ntrain, s_ref=symbols(n); fb_new=symbols(n); else, s_ref=d; fb_new=d; end
        e=s_ref-y; err(n)=e;
        Pu=P*u; kap=Pu/(lambda+u'*Pu);
        w=w+kap*conj(e); P=(P-kap*(Pu'))/lambda;
        fb=[fb_new; fb(1:end-1)];
    end
end
