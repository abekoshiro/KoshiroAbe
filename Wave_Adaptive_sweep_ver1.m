%% ============================================================
%  Wave_Adaptive_sweep_ver1.m
%  時変チャネル(波ドップラー)で【適応DFE(凍結OFF)＋λスイープ】を評価。
%
%  目的:
%   ver5の「学習後凍結」は静的最適で、波環境では時変OSNRが崩壊した（~4dB）。
%   チャネルがパケット内で老化（|c|~0.64/ステップ）するなら、学習後も
%   係数を動かし続ける【追従型RLS】が要る。忘却係数 λ を小さくするほど
%   速く追従するが、判定誤りに過剰反応しやすい。最適な λ を探す。
%
%  比較:
%   ・静的天井      : 波なし(snapshot1固定)での性能 … 到達可能な上限
%   ・時変＋凍結    : ver5そのまま(FREEZE=true) … 現状の崩壊ベースライン
%   ・時変＋適応λ   : FREEZE=false で λ を振る … どこまで取り戻せるか
%
%  効率: 重い処理(時変畳み込み/PTR/ビーム抽出)は試行ごとに1回。
%        同じビーム X を各等化器設定で等化し直すだけ。
%
%  依存: read_arrivals_asc.m, rrc_time_domain.m
% ============================================================
arr_base = 'JamTank_Arr';
n_snap   = 8;
f_wave   = 1.0;  T_wave = 1/f_wave;     % ★実験条件に合わせる

%% 信号・DFE（ver5準拠）
Rs=5000; Sps=16; Fs=Rs*Sps; numSymbols=2000;
alpha=0.25; Ncycle=max(20,ceil(4.5/alpha)); SNR_dB=20; Lfrac=16;
Nf=20; Nb=15; delta=0.01; mu_pll=0.01; Ntrain=600; Ntrials=10;

%% ★λスイープ範囲（適応RLS, FREEZE=false）
lam_list = [0.95 0.96 0.97 0.98 0.99 0.995 0.999];
lam_frozen = 0.999;   % 凍結・静的基準で使う学習時λ

%% 全スナップショット読込
Arr_all=cell(n_snap,1);
for k=1:n_snap
    fk=sprintf('%s_wave%02d.arr',arr_base,k);
    if ~isfile(fk), error('見つかりません: %s',fk); end
    [Arr_all{k},~]=read_arrivals_asc(fk);
end
Nrr=size(Arr_all{1},1); Nrz=size(Arr_all{1},2); Nsd=size(Arr_all{1},3); NUM_RX=Nrr*Nrz;
K_SUB=min(8,NUM_RX); ch_sub=min(K_SUB,ceil((1:NUM_RX)'/ceil(NUM_RX/K_SUB)))';
ch_rr=zeros(NUM_RX,1); ch_rz=zeros(NUM_RX,1); cc=0;
for irr=1:Nrr, for irz=1:Nrz, cc=cc+1; ch_rr(cc)=irr; ch_rz(cc)=irz; end, end
fprintf('ARR %dch 音源%d スナップ%d | f_wave=%.2fHz(0.4sパケット=%.2f周期) | λスイープ%d点\n\n',...
        NUM_RX,Nsd,n_snap,f_wave,(numSymbols/Rs)/T_wave,numel(lam_list));

%% RRC
rrc=rrc_time_domain(Rs,Ncycle,alpha,Fs); bad=~isfinite(rrc);
if any(bad), val=alpha/sqrt(2)*((1+2/pi)*sin(pi/(4*alpha))+(1-2/pi)*cos(pi/(4*alpha))); rrc(bad)=val; end
rrc=rrc(:)/sqrt(sum(abs(rrc).^2)); gd=Ncycle*Sps;

Nlam=numel(lam_list);
% OSNR格納: 各試行のdB値
O_static=nan(Nsd,Ntrials);            % 静的天井
O_frozen=nan(Nsd,Ntrials);            % 時変＋凍結
O_adapt =nan(Nsd,Ntrials,Nlam);       % 時変＋適応λ
syms_best=cell(Nsd,1); lam_best=nan(Nsd,1);

for idx_sd=1:Nsd
    %% 共通最小遅延・CIR長
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

    %% 受信機×スナップショットのCIR事前構築
    Hcir=cell(NUM_RX,1);
    for ch=1:NUM_RX
        Hk=zeros(Lh,n_snap);
        for k=1:n_snap
            amp=Arr_all{k}(ch_rr(ch),ch_rz(ch),idx_sd).A;
            dl =real(Arr_all{k}(ch_rr(ch),ch_rz(ch),idx_sd).delay);
            v=(dl>0); amp=amp(v); dl=dl(v); dl=dl-gmin; v=(dl>=0); amp=amp(v); dl=dl(v);
            if isempty(dl), continue; end
            Hk(:,k)=build_cir_frac(amp,dl*Fs+Lfrac,Lh,Lfrac);
        end
        Hcir{ch}=Hk;
    end
    Lblk=max(Sps,round(Fs*T_wave/64));
    ti=(Ntrain+1):numSymbols;

    for it=1:Ntrials
        bits=randi([0 1],2*numSymbols,1);
        symbols=(2*bits(1:2:end)-1)+1i*(2*bits(2:2:end)-1);
        up=zeros(numSymbols*Sps,1); up(1:Sps:end)=symbols; tx=conv(up,rrc); Ltx=length(tx);
        plen=Ltx+3*Lh+4*length(rrc);
        sub_ST=zeros(plen,K_SUB); sub_TV=zeros(plen,K_SUB); sub_ac=zeros(plen,K_SUB);

        for ch=1:NUM_RX
            Hk=Hcir{ch}; if all(Hk(:)==0), continue; end
            h_probe=Hk(:,1); g=conj(flipud(h_probe));
            rx_ST=fftconv(tx,h_probe);
            rx_TV=tv_channel(tx,Hk,Fs,T_wave,Lblk,n_snap);
            sp=mean(abs(rx_ST).^2); np=sp/10^(SNR_dB/10);
            L=max(length(rx_ST),length(rx_TV)); nz=sqrt(np/2)*(randn(L,1)+1i*randn(L,1));
            rx_ST(end+1:L)=0; rx_TV(end+1:L)=0; rx_ST=rx_ST+nz; rx_TV=rx_TV+nz;
            zST=fftconv(fftconv(rx_ST,rrc),g); zTV=fftconv(fftconv(rx_TV,rrc),g);
            hac=fftconv(h_probe,g); si=ch_sub(ch);
            LzS=min(length(zST),plen); sub_ST(1:LzS,si)=sub_ST(1:LzS,si)+zST(1:LzS);
            LzT=min(length(zTV),plen); sub_TV(1:LzT,si)=sub_TV(1:LzT,si)+zTV(1:LzT);
            La =min(length(hac),plen); sub_ac(1:La,si) =sub_ac(1:La,si) +hac(1:La);
        end

        %% ビーム抽出（静的・時変）— 試行ごとに1回
        X_ST=get_beams(sub_ST,sub_ac,ch_sub,K_SUB,Sps,numSymbols,symbols,gd);
        X_TV=get_beams(sub_TV,sub_ac,ch_sub,K_SUB,Sps,numSymbols,symbols,gd);

        %% 静的天井（凍結でOK）
        if ~isempty(X_ST)
            o=equalize_osnr(X_ST,symbols,Nf,Nb,lam_frozen,delta,Ntrain,mu_pll,true,ti);
            O_static(idx_sd,it)=o;
        end
        %% 時変＋凍結（崩壊ベースライン）
        if ~isempty(X_TV)
            O_frozen(idx_sd,it)=equalize_osnr(X_TV,symbols,Nf,Nb,lam_frozen,delta,Ntrain,mu_pll,true,ti);
            %% 時変＋適応λスイープ（同じXを等化し直すだけ）
            besto=-inf; bestc=[];
            for li=1:Nlam
                [o,c]=equalize_osnr(X_TV,symbols,Nf,Nb,lam_list(li),delta,Ntrain,mu_pll,false,ti);
                O_adapt(idx_sd,it,li)=o;
                if o>besto, besto=o; bestc=c; lam_best(idx_sd)=lam_list(li); end
            end
            if ~isempty(bestc), syms_best{idx_sd}=bestc; end
        end
    end
end

%% ==== 集計・出力 ====
disp('==========================================================================');
fprintf(' 適応DFE λスイープ（時変チャネル, f_wave=%.2fHz, %d試行）\n',f_wave,Ntrials);
disp('==========================================================================');
for idx_sd=1:Nsd
    oceil=median(O_static(idx_sd,:),'omitnan');
    ofrz =median(O_frozen(idx_sd,:),'omitnan');
    fprintf('【音源%d】 静的天井 %.1fdB | 時変＋凍結 %.1fdB（崩壊）\n',idx_sd,oceil,ofrz);
    fprintf('  時変＋適応λ:\n');
    oad=nan(Nlam,1);
    for li=1:Nlam
        oad(li)=median(O_adapt(idx_sd,:,li),'omitnan');
        fprintf('    λ=%.3f : %6.2f dB\n',lam_list(li),oad(li));
    end
    [bo,bi]=max(oad);
    fprintf('  → 最良 λ=%.3f で %.2f dB（凍結比 %+.2f dB, 天井まで残り %.2f dB）\n\n',...
            lam_list(bi),bo,bo-ofrz,oceil-bo);
end

%% ==== 図: OSNR vs λ ====
figure('Name','適応DFE λスイープ（時変チャネル）','Position',[100 100 460*Nsd 420]);
for idx_sd=1:Nsd
    oceil=median(O_static(idx_sd,:),'omitnan'); ofrz=median(O_frozen(idx_sd,:),'omitnan');
    oad=arrayfun(@(li) median(O_adapt(idx_sd,:,li),'omitnan'),1:Nlam);
    subplot(1,Nsd,idx_sd);
    plot(lam_list,oad,'-o','LineWidth',1.8,'MarkerFaceColor','w'); hold on;
    yline(oceil,'--','静的天井','Color',[0.1 0.5 0.1],'LineWidth',1.4,'LabelHorizontalAlignment','left');
    yline(ofrz,'--','時変＋凍結','Color',[0.85 0.3 0.1],'LineWidth',1.4,'LabelHorizontalAlignment','left');
    grid on; xlabel('忘却係数 λ（小=速い追従）'); ylabel('OSNR [dB]（中央値）');
    title(sprintf('音源%d: 適応DFEの回復',idx_sd));
    legend('時変＋適応λ','Location','best');
end
sgtitle(sprintf('波ドップラー下の等化器: 凍結は崩壊 → 追従で回復（f_{wave}=%.2fHz）',f_wave),'FontWeight','bold');

%% ================= ローカル関数 =================
function rx = tv_channel(tx, Hc, Fs, T_wave, Lblk, n_snap)
    Ltx=length(tx); Lh=size(Hc,1); rx=zeros(Ltx+Lh-1,1); nblk=ceil(Ltx/Lblk);
    for m=1:nblk
        a=(m-1)*Lblk+1; b=min(m*Lblk,Ltx); tc=((a+b)/2)/Fs;
        phi=2*pi*mod(tc,T_wave)/T_wave; h=cir_at_phase(Hc,phi,n_snap);
        c=fftconv(tx(a:b),h); rx(a:a+length(c)-1)=rx(a:a+length(c)-1)+c;
    end
end
function h = cir_at_phase(Hc, phi, n_snap)
    f=mod(phi,2*pi)/(2*pi)*n_snap; i0=floor(f); frac=f-i0;
    a=mod(i0,n_snap)+1; b=mod(i0+1,n_snap)+1; h=(1-frac)*Hc(:,a)+frac*Hc(:,b);
end
function X = get_beams(sub_ptr,sub_ac,ch_sub,K_SUB,Sps,numSymbols,symbols,gd)
    X=zeros(numSymbols,K_SUB); valid=false(K_SUB,1);
    for s=1:K_SUB
        if all(sub_ac(:,s)==0), continue; end
        [~,pk]=max(abs(sub_ac(:,s)));
        beam=find_sample(sub_ptr(:,s),pk,Sps,numSymbols,symbols,3*gd+4*Sps);
        if isempty(beam), continue; end
        a=(symbols'*beam)/(symbols'*symbols);
        if abs(a)<eps||any(~isfinite(beam)), continue; end
        if mean(abs(symbols-beam/a).^2)<3, X(:,s)=beam/a; valid(s)=true; end
    end
    X=X(:,valid);
end
function [osnr,comp] = equalize_osnr(X,symbols,Nf,Nb,lambda,delta,Ntrain,mu_pll,freeze,ti)
    out=mc_dfe(X,symbols,Nf,Nb,lambda,delta,Ntrain,mu_pll,freeze);
    if all(isfinite(out))
        [osnr,~,comp]=eval_region(symbols,out,ti);
    else
        osnr=0; comp=zeros(numel(symbols),1);
    end
end
function h = build_cir_frac(amp, cs, Lh, L)
    h=zeros(Lh,1);
    for p=1:numel(amp)
        c=cs(p); m=(floor(c)-L:floor(c)+L).'; ok=m>=0 & m<=Lh-1; m=m(ok); if isempty(m), continue; end
        x=m-c; w=0.5*(1+cos(pi*x/(L+1))); h(m+1)=h(m+1)+amp(p)*(sinc(x).*w);
    end
end
function y = fftconv(a, b), n=length(a)+length(b)-1; y=ifft(fft(a(:),n).*fft(b(:),n)); end
function sym = find_sample(sig, start, Sps, numSymbols, symbols, maxoff)
    bestd=0; best=-inf;
    for d=0:Sps:maxoff
        o=start+d; idx=o+(0:numSymbols-1)*Sps; if idx(1)<1||idx(end)>length(sig), continue; end
        v=osnr_of(sig(idx),symbols); if v>best, best=v; bestd=d; end
    end
    sym=[]; best=-inf;
    for d=bestd-Sps:bestd+Sps
        o=start+d; idx=o+(0:numSymbols-1)*Sps; if idx(1)<1||idx(end)>length(sig), continue; end
        cand=sig(idx); v=osnr_of(cand,symbols); if v>best, best=v; sym=cand; end
    end
end
function v = osnr_of(cand, symbols)
    g0=(symbols'*cand)/(symbols'*symbols); if abs(g0)<eps, v=-inf; return; end
    cc=cand/g0; v=sum(abs(symbols).^2)/max(sum(abs(symbols-cc).^2),eps);
end
function out = mc_dfe(X, symbols, Nf, Nb, lambda, delta, Ntrain, mu_pll, freeze)
    if nargin<9, freeze=true; end
    [N,K]=size(X); Ntap=K*Nf+Nb; w=zeros(Ntap,1);
    for k=1:K, w((k-1)*Nf+1)=1/K; end
    P=(1/delta)*eye(Ntap); ff=zeros(Nf,K); fb=zeros(Nb,1); theta=0; out=zeros(N,1);
    for n=1:N
        ff=[X(n,:); ff(1:end-1,:)]; ffvec=reshape(ff,[],1)*exp(-1i*theta);
        u=[ffvec; -fb]; y=w'*u; out(n)=y; d=sign(real(y))+1i*sign(imag(y));
        if n<=Ntrain, s_ref=symbols(n); fb_new=symbols(n); else, s_ref=d; fb_new=d; end
        if (n<=Ntrain) || ~freeze
            e=s_ref-y; Pu=P*u; kap=Pu/(lambda+u'*Pu); w=w+kap*conj(e); P=(P-kap*(Pu'))/lambda;
        end
        phi=imag(y*conj(s_ref)); theta=theta+mu_pll*phi; fb=[fb_new; fb(1:end-1)];
    end
end
function [osnr, ber, comp, psig, perr] = eval_region(symbols, raw, idx)
    s=symbols(idx); r=raw(idx); a=(s'*r)/(s'*s); c=r/a;
    psig=mean(abs(s).^2); perr=mean(abs(s-c).^2); osnr=10*log10(psig/perr);
    ber=mean([real(c)>0,imag(c)>0]~=[real(s)>0,imag(s)>0],'all'); comp=raw/a;
end
