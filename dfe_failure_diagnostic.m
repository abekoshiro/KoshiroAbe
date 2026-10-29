%% ============================================================
%  dfe_failure_diagnostic.m
%  RLSマルチチャネルDFEの「失敗試行」の収束を診断する。
%   ・ver3と同じRLSで複数試行を実行
%   ・各試行で MC-DFE の 誤差|e[n]| と 重みノルム‖w[n]‖ を記録
%   ・失敗試行(OSNR低)と成功試行を取り出し比較:
%       - MSE学習曲線 |e[n]|^2 : いつ崩れるか
%       - ‖w[n]‖ の推移      : RLS数値発散なら急増 → 専門家示唆の確認
%  依存: rrc_time_domain.m
% ============================================================
arr_file = 'JamTank_Arr_ver3.arr';    % ★失敗が起きる環境(199ch等)に合わせる
[Arr, Pos] = read_arrivals_asc(arr_file);
Nrr=size(Arr,1); Nrz=size(Arr,2); Nsd=size(Arr,3); NUM_RX=Nrr*Nrz;
fprintf('ARR: %dch, 音源%d\n', NUM_RX, Nsd);

%% ==== 設定（ver3のRLS条件を再現。失敗が出やすいα=0.5）====
idx_sd = 1;
Rs=5000; Sps=16; Fs=Rs*Sps; numSymbols=2000;
alpha=0.5; Ncycle=max(20,ceil(4.5/alpha)); SNR_dB=20; Lfrac=16;
K_SUB=min(8,NUM_RX); ch_sub=min(K_SUB,ceil((1:NUM_RX)'/ceil(NUM_RX/K_SUB)))';
Nf=20; Nb=15; lambda=0.999; delta=0.01; mu_pll=0.01; Ntrain=600;
Ntrials=40;                 % 失敗を捕まえるため多めに
fail_thr=10;                % OSNR<10dB を失敗とみなす

ch_rr=zeros(NUM_RX,1); ch_rz=zeros(NUM_RX,1); cc=0;
for irr=1:Nrr, for irz=1:Nrz, cc=cc+1; ch_rr(cc)=irr; ch_rz(cc)=irz; end, end

rrc=rrc_time_domain(Rs,Ncycle,alpha,Fs);
bad=~isfinite(rrc);
if any(bad), val=alpha/sqrt(2)*((1+2/pi)*sin(pi/(4*alpha))+(1-2/pi)*cos(pi/(4*alpha))); rrc(bad)=val; end
rrc=rrc(:)/sqrt(sum(abs(rrc).^2)); gd=Ncycle*Sps;

%% 共通最小遅延・CIR長
gmin=inf; Lmax=0;
for ch=1:NUM_RX, d=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay); d=d(d>0); if ~isempty(d), gmin=min(gmin,min(d)); end, end
for ch=1:NUM_RX, d=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay); d=d(d>0)-gmin; d=d(d>=0); if ~isempty(d), Lmax=max(Lmax,round(max(d)*Fs)+1); end, end
Lh=Lmax+2*Lfrac+4;

%% ==== 試行ループ ====
osnr_all=nan(Ntrials,1);
fail=struct('has',false); succ=struct('has',false);
nfail=0;
for it=1:Ntrials
    bits=randi([0 1],2*numSymbols,1);
    symbols=(2*bits(1:2:end)-1)+1i*(2*bits(2:2:end)-1);
    up=zeros(numSymbols*Sps,1); up(1:Sps:end)=symbols; tx=conv(up,rrc); Ltx=length(tx);
    ptr_buf_len=Ltx+3*Lh+4*length(rrc);
    sub_ptr=zeros(ptr_buf_len,K_SUB); sub_ac=zeros(ptr_buf_len,K_SUB);
    for ch=1:NUM_RX
        amp=Arr(ch_rr(ch),ch_rz(ch),idx_sd).A; dl=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay);
        v=(dl>0); amp=amp(v); dl=dl(v); dl=dl-gmin; v=(dl>=0); amp=amp(v); dl=dl(v);
        if isempty(dl), continue; end
        cs=dl*Fs+Lfrac; h=build_cir_frac(amp,cs,Lh,Lfrac);
        rxBase=fftconv(tx,h); sp=mean(abs(rxBase).^2); np=sp/10^(SNR_dB/10);
        rxn=rxBase+sqrt(np/2)*(randn(size(rxBase))+1i*randn(size(rxBase)));
        rxmf=fftconv(rxn,rrc); g=conj(flipud(h));
        z=fftconv(rxmf,g); hac=fftconv(h,g);
        Lz=min(length(z),ptr_buf_len); La=min(length(hac),ptr_buf_len);
        si=ch_sub(ch); sub_ptr(1:Lz,si)=sub_ptr(1:Lz,si)+z(1:Lz); sub_ac(1:La,si)=sub_ac(1:La,si)+hac(1:La);
    end
    X=zeros(numSymbols,K_SUB); valid=false(K_SUB,1);
    for s=1:K_SUB
        if all(sub_ac(:,s)==0), continue; end
        [~,pk]=max(abs(sub_ac(:,s)));
        beam=find_sample(sub_ptr(:,s),pk,Sps,numSymbols,symbols,3*gd+4*Sps);
        if isempty(beam), continue; end
        a=(symbols'*beam)/(symbols'*symbols);
        if abs(a)<eps||any(~isfinite(beam)), continue; end
        eb=mean(abs(symbols-beam/a).^2);
        if eb<3, X(:,s)=beam/a; valid(s)=true; end
    end
    X=X(:,valid); if size(X,2)==0, osnr_all(it)=-30; nfail=nfail+1; continue; end
    ti=(Ntrain+1):numSymbols;

    %% RLS MC-DFE（誤差・重みノルムを記録）
    [mc_out,errh,wnormh]=mc_dfe_log(X,symbols,Nf,Nb,lambda,delta,Ntrain,mu_pll);
    if all(isfinite(mc_out))
        [osnr,~,comp]=eval_region(symbols,mc_out,ti);
    else
        osnr=-30; comp=zeros(numSymbols,1);
    end
    osnr_all(it)=osnr;

    isfail = (osnr<fail_thr) || ~all(isfinite(mc_out));
    if isfail, nfail=nfail+1; end
    if isfail && ~fail.has
        fail.has=true; fail.err=errh; fail.wn=wnormh; fail.comp=comp; fail.osnr=osnr; fail.it=it;
    elseif ~isfail && ~succ.has
        succ.has=true; succ.err=errh; succ.wn=wnormh; succ.comp=comp; succ.osnr=osnr; succ.it=it;
    end
end

fprintf('\n=== 結果（%d試行, α=%.2f）===\n', Ntrials, alpha);
fprintf('  失敗率: %d/%d = %.0f%%（OSNR<%ddB）\n', nfail, Ntrials, 100*nfail/Ntrials, fail_thr);
fprintf('  OSNR中央値: %.2f dB\n', median(osnr_all,'omitnan'));
if fail.has, fprintf('  失敗試行#%d: OSNR=%.1fdB, ‖w‖最大=%.3e\n', fail.it, fail.osnr, max(fail.wn)); end
if succ.has, fprintf('  成功試行#%d: OSNR=%.1fdB, ‖w‖最大=%.3e\n', succ.it, succ.osnr, max(succ.wn)); end
if fail.has && succ.has
    fprintf('  → ‖w‖最大の比(失敗/成功)= %.1f倍', max(fail.wn)/max(succ.wn));
    if max(fail.wn) > 10*max(succ.wn)
        fprintf('　★失敗時に重みが発散＝RLS数値不安定が原因（専門家示唆を確認）\n');
    else
        fprintf('　重み発散は顕著でない＝誤り伝播/収束不足の可能性\n');
    end
end

%% ==== 描画 ====
figure('Name','RLS失敗試行の収束診断','Position',[80 80 1100 420]);
sm=@(x) movmean(x,31);

subplot(1,3,1); hold on; grid on;
if succ.has, plot(1:numSymbols,10*log10(sm(succ.err.^2)+eps),'b','LineWidth',1.4,'DisplayName',sprintf('成功 %.1fdB',succ.osnr)); end
if fail.has, plot(1:numSymbols,10*log10(sm(fail.err.^2)+eps),'r','LineWidth',1.4,'DisplayName',sprintf('失敗 %.1fdB',fail.osnr)); end
xline(Ntrain,'k--','訓練終了'); xlabel('シンボル番号 n'); ylabel('MSE=|e|^2 [dB]');
title('MSE学習曲線'); legend('Location','best');

subplot(1,3,2); hold on; grid on;
if succ.has, semilogy(1:numSymbols,succ.wn,'b','LineWidth',1.4,'DisplayName','成功'); end
if fail.has, semilogy(1:numSymbols,fail.wn,'r','LineWidth',1.4,'DisplayName','失敗'); end
set(gca,'YScale','log'); xline(Ntrain,'k--'); xlabel('シンボル番号 n'); ylabel('重みノルム ‖w‖');
title('重みノルムの推移（発散の有無）'); legend('Location','best');

subplot(1,3,3); hold on; grid on; axis equal;
qpsk=[1+1i,1-1i,-1+1i,-1-1i];
if fail.has, plot(real(fail.comp((Ntrain+1):end)),imag(fail.comp((Ntrain+1):end)),'.','Color',[0.9 0.4 0.4],'MarkerSize',3); end
if succ.has, plot(real(succ.comp((Ntrain+1):end)),imag(succ.comp((Ntrain+1):end)),'.','Color',[0.4 0.6 1],'MarkerSize',3); end
plot(real(qpsk),imag(qpsk),'r+','MarkerSize',12,'LineWidth',2);
L=5; xlim([-L L]); ylim([-L L]); xlabel('I'); ylabel('Q'); title('星座図（青:成功 赤:失敗）');
sgtitle(sprintf('RLSマルチチャネルDFE 失敗試行の収束診断（α=%.2f, %dch）',alpha,NUM_RX),'FontWeight','bold');

%% ================= ローカル関数 =================
function [out,errh,wnormh] = mc_dfe_log(X, symbols, Nf, Nb, lambda, delta, Ntrain, mu_pll)
    [N,K]=size(X); Ntap=K*Nf+Nb; w=zeros(Ntap,1);
    for k=1:K, w((k-1)*Nf+1)=1/K; end
    P=(1/delta)*eye(Ntap); ff=zeros(Nf,K); fb=zeros(Nb,1); theta=0;
    out=zeros(N,1); errh=zeros(N,1); wnormh=zeros(N,1);
    for n=1:N
        ff=[X(n,:); ff(1:end-1,:)];
        ffvec=reshape(ff,[],1)*exp(-1i*theta);
        u=[ffvec; -fb]; y=w'*u; out(n)=y;
        d=sign(real(y))+1i*sign(imag(y));
        if n<=Ntrain, s_ref=symbols(n); fb_new=symbols(n); else, s_ref=d; fb_new=d; end
        e=s_ref-y; errh(n)=abs(e);
        Pu=P*u; kap=Pu/(lambda+u'*Pu);
        w=w+kap*conj(e); P=(P-kap*(Pu'))/lambda;
        wnormh(n)=norm(w);
        phi=imag(y*conj(s_ref)); theta=theta+mu_pll*phi;
        fb=[fb_new; fb(1:end-1)];
    end
end

function h = build_cir_frac(amp, cs, Lh, L)
    h=zeros(Lh,1);
    for p=1:numel(amp)
        c=cs(p); m=(floor(c)-L:floor(c)+L).';
        ok=m>=0 & m<=Lh-1; m=m(ok); if isempty(m), continue; end
        x=m-c; w=0.5*(1+cos(pi*x/(L+1))); h(m+1)=h(m+1)+amp(p)*(sinc(x).*w);
    end
end
function y = fftconv(a, b), n=length(a)+length(b)-1; y=ifft(fft(a(:),n).*fft(b(:),n)); end
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
function [osnr, ber, comp] = eval_region(symbols, raw, idx)
    s=symbols(idx); r=raw(idx); a=(s'*r)/(s'*s); c=r/a;
    osnr=10*log10(mean(abs(s).^2)/mean(abs(s-c).^2));
    ber=mean([real(c)>0,imag(c)>0]~=[real(s)>0,imag(s)>0],'all'); comp=raw/a;
end
