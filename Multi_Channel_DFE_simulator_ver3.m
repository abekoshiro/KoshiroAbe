%% ============================================================
%  Multi_Channel_DFE_simulator_ver3.m
%  ver2（RRC導入）に【遅延の丸め誤差の解消：sinc分数遅延】を追加した版。
%
%  専門家指摘①への対応（時間領域のまま）:
%   ・各パスを round(τ*Fs) で最寄りサンプルにスナップせず、
%     窓付きsinc補間で「厳密な小数サンプル位置」に配置（分数遅延）
%   ・分数遅延CIR h を受信生成(rx=tx⊛h)とPTR整合(g=conj(flip h))で共用 → 整合維持
%   ・正確な遅延・位相関係を保持 → PTR集束の精度向上
%
%  他はver2を踏襲: 送信RRC整形 / 受信整合RRC / 帯域内SNR / PTR焦点整列(共通長)
%  高速化: 大きな畳み込みはFFT (fftconv)
%  依存: rrc_time_domain.m
% ============================================================
arr_file = 'JamTank_Arr_ver3.arr';    % ★使用中のARRに合わせる
[Arr, Pos] = read_arrivals_asc(arr_file);
Nrr=size(Arr,1); Nrz=size(Arr,2); Nsd=size(Arr,3); NUM_RX=Nrr*Nrz;
fprintf('ARR: 距離%d × 深度%d = %dch, 音源%d\n', Nrr,Nrz,NUM_RX,Nsd);

%% ==== 信号パラメータ ====
Rs=5000; Sps=16; Fs=Rs*Sps; numSymbols=2000;   % BW=Rs, Fs=80kHz
alpha  = 0.25;    % ロールオフ率
Ncycle = 20;      % RRCスパン
SNR_dB = 20;      % 帯域内SNR
Lfrac  = 16;      % ★分数遅延sincの半長（大きいほど正確, 重い）

%% ==== サブアレイ・DFE ====
K_SUB = min(8, NUM_RX);
ch_sub = min(K_SUB, ceil((1:NUM_RX)'/ceil(NUM_RX/K_SUB)))';
Nf=20; Nb=15; lambda=0.999; delta=0.01; mu_pll=0.01; Ntrain=600;
Ntrials = 10;

fprintf('RRC α=%.2f Ncycle=%d | 分数遅延Lfrac=%d | K=%d Nf=%d/ch Nb=%d | Rs=%d %d試行\n\n', ...
        alpha,Ncycle,Lfrac,K_SUB,Nf,Nb,Rs,Ntrials);

ch_rr=zeros(NUM_RX,1); ch_rz=zeros(NUM_RX,1); cc=0;
for irr=1:Nrr, for irz=1:Nrz, cc=cc+1; ch_rr(cc)=irr; ch_rz(cc)=irz; end, end

%% ==== RRC（送受共用）====
rrc = rrc_time_domain(Rs, Ncycle, alpha, Fs);
bad=~isfinite(rrc);
if any(bad)
    val=alpha/sqrt(2)*((1+2/pi)*sin(pi/(4*alpha))+(1-2/pi)*cos(pi/(4*alpha)));
    rrc(bad)=val; fprintf('  （RRC特異点 %d 点を補修）\n',nnz(bad));
end
rrc=rrc(:)/sqrt(sum(abs(rrc).^2)); gd=Ncycle*Sps;

SCM_Ps=zeros(Nsd,Ntrials); SCM_Pe=zeros(Nsd,Ntrials); SCM_B=zeros(Nsd,Ntrials);
MC_Ps =zeros(Nsd,Ntrials); MC_Pe =zeros(Nsd,Ntrials); MC_B =zeros(Nsd,Ntrials);
SCM_syms=cell(Nsd,1); MC_syms=cell(Nsd,1);

for idx_sd=1:Nsd
    %% 共通最小遅延・共通CIR長（分数遅延の裾を考慮して Lfrac 分の余裕）
    gmin=inf; Lmax=0;
    for ch=1:NUM_RX
        d=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay); d=d(d>0);
        if ~isempty(d), gmin=min(gmin,min(d)); end
    end
    for ch=1:NUM_RX
        d=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay); d=d(d>0)-gmin; d=d(d>=0);
        if ~isempty(d), Lmax=max(Lmax, round(max(d)*Fs)+1); end
    end
    Lh = Lmax + 2*Lfrac + 4;      % 共通CIR長

    for it=1:Ntrials
        %% 送信 QPSK → RRC
        bits=randi([0 1],2*numSymbols,1);
        symbols=(2*bits(1:2:end)-1)+1i*(2*bits(2:2:end)-1);
        up=zeros(numSymbols*Sps,1); up(1:Sps:end)=symbols;
        tx=conv(up,rrc); Ltx=length(tx);

        ptr_buf_len=Ltx+3*Lh+4*length(rrc);
        sub_ptr=zeros(ptr_buf_len,K_SUB); sub_ac=zeros(ptr_buf_len,K_SUB);

        for ch=1:NUM_RX
            amp=Arr(ch_rr(ch),ch_rz(ch),idx_sd).A; dl=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay);
            v=(dl>0); amp=amp(v); dl=dl(v); dl=dl-gmin; v=(dl>=0); amp=amp(v); dl=dl(v);
            if isempty(dl), continue; end

            %% ★分数遅延CIR（sinc補間, +Lfracで先頭パスの裾切れ防止）
            cs = dl*Fs + Lfrac;
            h  = build_cir_frac(amp, cs, Lh, Lfrac);

            %% チャネル適用（rx=tx⊛h, 厳密な小数遅延）→ 帯域内SNR → 受信整合RRC
            rxBase=fftconv(tx,h);
            sp=mean(abs(rxBase).^2); np=sp/10^(SNR_dB/10);
            rxn=rxBase+sqrt(np/2)*(randn(size(rxBase))+1i*randn(size(rxBase)));
            rxmf=fftconv(rxn,rrc);

            %% PTR（同じ分数遅延CIRで整合フィルタ）
            g=conj(flipud(h));
            z=fftconv(rxmf,g); hac=fftconv(h,g);
            Lz=min(length(z),ptr_buf_len); La=min(length(hac),ptr_buf_len);
            si=ch_sub(ch);
            sub_ptr(1:Lz,si)=sub_ptr(1:Lz,si)+z(1:Lz);
            sub_ac(1:La,si) =sub_ac(1:La,si)+hac(1:La);
        end

        %% ビーム抽出（焦点を広域探索）＋正規化
        X=zeros(numSymbols,K_SUB); valid=false(K_SUB,1); err_beam=zeros(K_SUB,1);
        for s=1:K_SUB
            if all(sub_ac(:,s)==0), continue; end
            [~,pk]=max(abs(sub_ac(:,s)));
            beam=find_sample(sub_ptr(:,s),pk,Sps,numSymbols,symbols,3*gd+4*Sps);
            if isempty(beam), continue; end
            a=(symbols'*beam)/(symbols'*symbols);
            if abs(a)<eps||any(~isfinite(beam)), continue; end
            X(:,s)=beam/a; valid(s)=true; err_beam(s)=mean(abs(symbols-X(:,s)).^2);
        end
        X=X(:,valid); err_beam=err_beam(valid); Kv=size(X,2);
        if Kv==0, continue; end
        ti=(Ntrain+1):numSymbols;

        %% 参考: SCM-MRC → 単一DFE
        scm=zeros(numSymbols,1); wsum=0;
        for k=1:Kv, w=1/err_beam(k); scm=scm+w*X(:,k); wsum=wsum+w; end
        scm=scm/wsum;
        scm_dfe=rls_dfe(scm,symbols,Nf,Nb,lambda,delta,Ntrain);
        if all(isfinite(scm_dfe))
            [~,SCM_B(idx_sd,it),cS,SCM_Ps(idx_sd,it),SCM_Pe(idx_sd,it)]=eval_region(symbols,scm_dfe,ti);
        else
            SCM_Ps(idx_sd,it)=1; SCM_Pe(idx_sd,it)=1; SCM_B(idx_sd,it)=0.5; cS=zeros(numSymbols,1);
        end

        %% 本命: マルチチャネルDFE
        mc_out=mc_dfe(X,symbols,Nf,Nb,lambda,delta,Ntrain,mu_pll);
        [~,MC_B(idx_sd,it),cM,MC_Ps(idx_sd,it),MC_Pe(idx_sd,it)]=eval_region(symbols,mc_out,ti);

        if it==Ntrials, SCM_syms{idx_sd}=cS; MC_syms{idx_sd}=cM; end
    end
end

%% ==== 集計・出力 ====
linO=@(Ps,Pe) 10*log10(mean(Ps,2)./mean(Pe,2)); mB=@(B) mean(B,2);
disp('==========================================================================');
fprintf(' 分数遅延版: SCM+単一DFE vs マルチチャネルDFE（%d試行, Rs=%d, α=%.2f, Lfrac=%d）\n', Ntrials,Rs,alpha,Lfrac);
disp('==========================================================================');
for idx_sd=1:Nsd
    oS=linO(SCM_Ps(idx_sd,:),SCM_Pe(idx_sd,:)); oM=linO(MC_Ps(idx_sd,:),MC_Pe(idx_sd,:));
    fprintf('【音源%d】\n',idx_sd);
    fprintf('  SCM-MRC + 単一DFE : OSNR %6.2f dB   BER %.4f\n', oS, mB(SCM_B(idx_sd,:)));
    fprintf('  マルチチャネルDFE : OSNR %6.2f dB   BER %.4f   （改善 %+.2f dB）\n\n', oM, mB(MC_B(idx_sd,:)), oM-oS);
end

%% ==== 星座図 ====
qpsk=[1+1i,1-1i,-1+1i,-1-1i];
figure('Name','分数遅延版: SCM vs MC-DFE','Position',[100 100 300*Nsd 560]);
for idx_sd=1:Nsd
    oS=linO(SCM_Ps(idx_sd,:),SCM_Pe(idx_sd,:)); oM=linO(MC_Ps(idx_sd,:),MC_Pe(idx_sd,:));
    subplot(2,Nsd,idx_sd);
    plot(real(SCM_syms{idx_sd}),imag(SCM_syms{idx_sd}),'.','Color',[0.1 0.6 0.9],'MarkerSize',3); hold on;
    plot(real(qpsk),imag(qpsk),'r+','MarkerSize',12,'LineWidth',2); axis equal; grid on;
    lim=max(3,ceil(max(abs([real(SCM_syms{idx_sd});imag(SCM_syms{idx_sd})])))); xlim([-lim lim]); ylim([-lim lim]);
    title(sprintf('SCM+単一DFE 音源%d\nOSNR=%.1fdB',idx_sd,oS));
    subplot(2,Nsd,Nsd+idx_sd);
    plot(real(MC_syms{idx_sd}),imag(MC_syms{idx_sd}),'.','Color',[0.85 0.3 0.1],'MarkerSize',3); hold on;
    plot(real(qpsk),imag(qpsk),'r+','MarkerSize',12,'LineWidth',2); axis equal; grid on;
    lim=max(3,ceil(max(abs([real(MC_syms{idx_sd});imag(MC_syms{idx_sd})])))); xlim([-lim lim]); ylim([-lim lim]);
    title(sprintf('マルチチャネルDFE 音源%d\nOSNR=%.1fdB',idx_sd,oM));
end
sgtitle(sprintf('分数遅延版（α=%.2f, Lfrac=%d, K=%d, Nf=%d/ch, Nb=%d, Rs=%d）',alpha,Lfrac,K_SUB,Nf,Nb,Rs),'FontWeight','bold');

%% ================= ローカル関数 =================
function h = build_cir_frac(amp, cs, Lh, L)
% 分数遅延CIR：各パスを窓付きsincで小数サンプル位置 cs(p) に配置
    h=zeros(Lh,1);
    for p=1:numel(amp)
        c=cs(p); m=(floor(c)-L:floor(c)+L).';
        ok=m>=0 & m<=Lh-1; m=m(ok); if isempty(m), continue; end
        x=m-c; w=0.5*(1+cos(pi*x/(L+1)));     % Hann窓
        h(m+1)=h(m+1)+amp(p)*(sinc(x).*w);    % sinc補間（分数遅延）
    end
end

function y = fftconv(a, b)
% FFTによる線形畳み込み（複素対応, 高速）
    n=length(a)+length(b)-1; y=ifft(fft(a(:),n).*fft(b(:),n));
end

function sym = find_sample(sig, start, Sps, numSymbols, symbols, maxoff)
    bestd=0; best=-inf;
    for d=0:Sps:maxoff
        o=start+d; idx=o+(0:numSymbols-1)*Sps;
        if idx(1)<1||idx(end)>length(sig), continue; end
        v=osnr_of(sig(idx),symbols); if v>best, best=v; bestd=d; end
    end
    sym=[]; best=-inf;
    for d=bestd-Sps:bestd+Sps
        o=start+d; idx=o+(0:numSymbols-1)*Sps;
        if idx(1)<1||idx(end)>length(sig), continue; end
        cand=sig(idx); v=osnr_of(cand,symbols); if v>best, best=v; sym=cand; end
    end
end

function v = osnr_of(cand, symbols)
    g0=(symbols'*cand)/(symbols'*symbols);
    if abs(g0)<eps, v=-inf; return; end
    cc=cand/g0; v=sum(abs(symbols).^2)/max(sum(abs(symbols-cc).^2),eps);
end

function out = mc_dfe(X, symbols, Nf, Nb, lambda, delta, Ntrain, mu_pll)
    [N,K]=size(X); Ntap=K*Nf+Nb; w=zeros(Ntap,1);
    for k=1:K, w((k-1)*Nf+1)=1/K; end
    P=(1/delta)*eye(Ntap); ff=zeros(Nf,K); fb=zeros(Nb,1); theta=0; out=zeros(N,1);
    for n=1:N
        ff=[X(n,:); ff(1:end-1,:)];
        ffvec=reshape(ff,[],1)*exp(-1i*theta);
        u=[ffvec; -fb]; y=w'*u; out(n)=y;
        d=sign(real(y))+1i*sign(imag(y));
        if n<=Ntrain, s_ref=symbols(n); fb_new=symbols(n); else, s_ref=d; fb_new=d; end
        e=s_ref-y; Pu=P*u; kap=Pu/(lambda+u'*Pu);
        w=w+kap*conj(e); P=(P-kap*(Pu'))/lambda;
        phi=imag(y*conj(s_ref)); theta=theta+mu_pll*phi;
        fb=[fb_new; fb(1:end-1)];
    end
end

function out = rls_dfe(x, symbols, Nf, Nb, lambda, delta, Ntrain)
    N=length(x); Ntap=Nf+Nb; w=zeros(Ntap,1); w(1)=1; P=(1/delta)*eye(Ntap);
    ff=zeros(Nf,1); fb=zeros(Nb,1); out=zeros(N,1);
    for n=1:N
        ff=[x(n); ff(1:end-1)]; u=[ff; -fb]; y=w'*u; out(n)=y;
        d=sign(real(y))+1i*sign(imag(y));
        if n<=Ntrain, s_ref=symbols(n); fb_new=symbols(n); else, s_ref=d; fb_new=d; end
        e=s_ref-y; Pu=P*u; kap=Pu/(lambda+u'*Pu);
        w=w+kap*conj(e); P=(P-kap*(Pu'))/lambda; fb=[fb_new; fb(1:end-1)];
    end
end

function [osnr, ber, comp, psig, perr] = eval_region(symbols, raw, idx)
    s=symbols(idx); r=raw(idx); a=(s'*r)/(s'*s); c=r/a;
    psig=mean(abs(s).^2); perr=mean(abs(s-c).^2); osnr=10*log10(psig/perr);
    ber=mean([real(c)>0,imag(c)>0]~=[real(s)>0,imag(s)>0],'all'); comp=raw/a;
end
