%% ============================================================
%  Wave_TimeVarying_simulator_ver1.m
%  【Level2】波で動く水面による「時変チャネル」でPTR+MC-DFEを評価する版。
%
%  考え方（凍結面スナップショット → 時変チャネル）:
%   ・run_wave_bellhop.m が作った 8枚の ARR（1波周期を8位相で刻んだもの）を、
%     実時間軸 t_k = (k-1)/n_snap * T_wave に並べる。
%   ・各受信機ごとに 8枚ぶんの分数遅延CIR h_k(τ) を作り、
%     送信時間に応じて波位相 φ(t)=2π·mod(t,T_wave)/T_wave で【CIRを補間】。
%   ・パケットをブロックに分け、各ブロックに「その時刻の補間CIR」を適用する
%     ブロック定常近似＋オーバーラップ加算で、時変畳み込みを実装。
%     → 水面が動くことによるドップラー／パケット内フェージングが再現される。
%
%  PTRの現実的モデル:
%   ・整合フィルタ g はパケット冒頭のチャネル推定（スナップショット1）で固定。
%   ・データ区間ではチャネルが φ(t) に沿ってドリフト → g とズレる
%     （＝ドップラーによるPTR集束劣化）。残差を MC-DFE(ver5凍結) で回復。
%
%  比較: 「静的（snapshot1固定）」vs「時変（ドップラーあり）」で OSNR を並べ、
%        波によるPTR劣化量と、MC-DFEがどこまで取り戻せるかを見る。
%
%  依存: read_arrivals_asc.m, rrc_time_domain.m
%  ★要設定: f_wave（造波機の周波数）。これでドップラーの大きさが決まる。
% ============================================================
arr_base = 'JamTank_Arr';
n_snap   = 8;                 % スナップショット数（make_wave_ati の n_snap と一致）

%% ==== ★波パラメータ（実験条件に合わせて必ず設定）====
f_wave = 1.0;                 % ★造波機の周波数 [Hz]（T_wave=1/f_wave）
                              %   不明なら分散関係から推定も可:
                              %   浅海 c=sqrt(g*depth), 深海 c=sqrt(g*lambda/2pi), T=lambda/c
T_wave = 1/f_wave;            % 波周期 [s]

%% ==== 信号パラメータ（ver5に合わせる）====
Rs=5000; Sps=16; Fs=Rs*Sps; numSymbols=2000;
alpha  = 0.25;
Ncycle = 20; Ncycle = max(Ncycle, ceil(4.5/alpha));
SNR_dB = 20;
Lfrac  = 16;

%% ==== サブアレイ・DFE（ver5準拠）====
Nf=20; Nb=15; lambda=0.999; delta=0.01; mu_pll=0.01; Ntrain=600;
Ntrials = 10; FREEZE = true;

%% ==== 全スナップショットの ARR を読込 ====
Arr_all = cell(n_snap,1);
for k=1:n_snap
    fk = sprintf('%s_wave%02d.arr', arr_base, k);
    if ~isfile(fk), error('見つかりません: %s （先に run_wave_bellhop.m）', fk); end
    [Arr_all{k}, Pos] = read_arrivals_asc(fk);
end
Nrr=size(Arr_all{1},1); Nrz=size(Arr_all{1},2); Nsd=size(Arr_all{1},3);
NUM_RX=Nrr*Nrz;
fprintf('ARR: 距離%d × 深度%d = %dch, 音源%d, スナップ%d枚\n', Nrr,Nrz,NUM_RX,Nsd,n_snap);
fprintf('波: f_wave=%.2fHz T_wave=%.3fs | パケット長=%.3fs (%.2f周期分)\n', ...
        f_wave, T_wave, numSymbols/Rs, (numSymbols/Rs)/T_wave);

K_SUB = min(8, NUM_RX);
ch_sub = min(K_SUB, ceil((1:NUM_RX)'/ceil(NUM_RX/K_SUB)))';
fprintf('K=%d Nf=%d/ch Nb=%d | α=%.2f Lfrac=%d | %d試行 FREEZE=%d\n\n', ...
        K_SUB,Nf,Nb,alpha,Lfrac,Ntrials,FREEZE);

ch_rr=zeros(NUM_RX,1); ch_rz=zeros(NUM_RX,1); cc=0;
for irr=1:Nrr, for irz=1:Nrz, cc=cc+1; ch_rr(cc)=irr; ch_rz(cc)=irz; end, end

%% ==== RRC（送受共用）====
rrc = rrc_time_domain(Rs, Ncycle, alpha, Fs);
bad=~isfinite(rrc);
if any(bad)
    val=alpha/sqrt(2)*((1+2/pi)*sin(pi/(4*alpha))+(1-2/pi)*cos(pi/(4*alpha)));
    rrc(bad)=val;
end
rrc=rrc(:)/sqrt(sum(abs(rrc).^2)); gd=Ncycle*Sps;

% 結果格納（stat=静的, tv=時変）
ST_Ps=zeros(Nsd,Ntrials); ST_Pe=zeros(Nsd,Ntrials);
TV_Ps=zeros(Nsd,Ntrials); TV_Pe=zeros(Nsd,Ntrials);
ST_B =zeros(Nsd,Ntrials); TV_B =zeros(Nsd,Ntrials);
ST_syms=cell(Nsd,1); TV_syms=cell(Nsd,1);

for idx_sd=1:Nsd
    %% 共通最小遅延・共通CIR長（全スナップショットで共通にそろえる）
    gmin=inf; Lmax=0;
    for k=1:n_snap
        for ch=1:NUM_RX
            d=real(Arr_all{k}(ch_rr(ch),ch_rz(ch),idx_sd).delay); d=d(d>0);
            if ~isempty(d), gmin=min(gmin,min(d)); end
        end
    end
    for k=1:n_snap
        for ch=1:NUM_RX
            d=real(Arr_all{k}(ch_rr(ch),ch_rz(ch),idx_sd).delay); d=d(d>0)-gmin; d=d(d>=0);
            if ~isempty(d), Lmax=max(Lmax, round(max(d)*Fs)+1); end
        end
    end
    Lh = Lmax + 2*Lfrac + 4;

    %% 受信機ごと・スナップショットごとの分数遅延CIR を事前構築: Hcir{ch} = Lh x n_snap
    Hcir = cell(NUM_RX,1);
    for ch=1:NUM_RX
        Hk = zeros(Lh,n_snap);
        for k=1:n_snap
            amp=Arr_all{k}(ch_rr(ch),ch_rz(ch),idx_sd).A;
            dl =real(Arr_all{k}(ch_rr(ch),ch_rz(ch),idx_sd).delay);
            v=(dl>0); amp=amp(v); dl=dl(v); dl=dl-gmin; v=(dl>=0); amp=amp(v); dl=dl(v);
            if isempty(dl), continue; end
            cs = dl*Fs + Lfrac;
            Hk(:,k) = build_cir_frac(amp, cs, Lh, Lfrac);
        end
        Hcir{ch}=Hk;
    end

    % ブロック長：波周期を細かく刻む（チャネル変化がブロック内で小さくなるよう）
    Lblk = max(Sps, round(Fs*T_wave/64));

    for it=1:Ntrials
        %% 送信 QPSK → RRC
        bits=randi([0 1],2*numSymbols,1);
        symbols=(2*bits(1:2:end)-1)+1i*(2*bits(2:2:end)-1);
        up=zeros(numSymbols*Sps,1); up(1:Sps:end)=symbols;
        tx=conv(up,rrc); Ltx=length(tx);

        ptr_buf_len=Ltx+3*Lh+4*length(rrc);
        sub_ST=zeros(ptr_buf_len,K_SUB); sub_TV=zeros(ptr_buf_len,K_SUB);
        sub_ac=zeros(ptr_buf_len,K_SUB);   % 焦点検出用（プローブCIRの自己相関）

        for ch=1:NUM_RX
            Hk=Hcir{ch}; if all(Hk(:)==0), continue; end
            h_probe=Hk(:,1);                 % プローブ=冒頭(snapshot1)のCIR
            g=conj(flipud(h_probe));          % PTR整合フィルタ（固定）

            %% 静的チャネル（snapshot1固定, ドップラーなし）
            rx_ST=fftconv(tx,h_probe);
            %% 時変チャネル（波でドリフト, ドップラーあり）
            rx_TV=tv_channel(tx,Hk,Fs,T_wave,Lblk,n_snap);

            %% 共通SNRで雑音（静的の受信電力を基準に同じ雑音レベル）
            sp=mean(abs(rx_ST).^2); np=sp/10^(SNR_dB/10);
            L=max(length(rx_ST),length(rx_TV));
            nz=sqrt(np/2)*(randn(L,1)+1i*randn(L,1));
            rx_ST(end+1:L)=0; rx_TV(end+1:L)=0;
            rx_ST=rx_ST+nz; rx_TV=rx_TV+nz;

            %% 受信整合RRC → PTR（同じ固定 g）
            zST=fftconv(fftconv(rx_ST,rrc),g);
            zTV=fftconv(fftconv(rx_TV,rrc),g);
            hac=fftconv(h_probe,g);

            si=ch_sub(ch);
            LzS=min(length(zST),ptr_buf_len); sub_ST(1:LzS,si)=sub_ST(1:LzS,si)+zST(1:LzS);
            LzT=min(length(zTV),ptr_buf_len); sub_TV(1:LzT,si)=sub_TV(1:LzT,si)+zTV(1:LzT);
            La =min(length(hac),ptr_buf_len); sub_ac(1:La,si) =sub_ac(1:La,si) +hac(1:La);
        end

        %% ビーム抽出（静的・時変で共通の焦点位置を使う）＋評価
        ti=(Ntrain+1):numSymbols;
        [cS,psS,peS,bS]=process_beams(sub_ST,sub_ac,ch_sub,K_SUB,Sps,numSymbols,symbols, ...
                                      gd,Nf,Nb,lambda,delta,Ntrain,mu_pll,FREEZE,ti);
        [cT,psT,peT,bT]=process_beams(sub_TV,sub_ac,ch_sub,K_SUB,Sps,numSymbols,symbols, ...
                                      gd,Nf,Nb,lambda,delta,Ntrain,mu_pll,FREEZE,ti);
        ST_Ps(idx_sd,it)=psS; ST_Pe(idx_sd,it)=peS; ST_B(idx_sd,it)=bS;
        TV_Ps(idx_sd,it)=psT; TV_Pe(idx_sd,it)=peT; TV_B(idx_sd,it)=bT;
        if any(isfinite(cS)), ST_syms{idx_sd}=cS; end
        if any(isfinite(cT)), TV_syms{idx_sd}=cT; end
    end
end

%% ==== 集計・出力 ====
linO=@(Ps,Pe) 10*log10(mean(Ps,2)./mean(Pe,2)); mB=@(B) mean(B,2);
disp('==========================================================================');
fprintf(' 時変チャネル版(Level2): 静的 vs 時変(ドップラーあり)  f_wave=%.2fHz\n', f_wave);
disp('==========================================================================');
for idx_sd=1:Nsd
    oS=linO(ST_Ps(idx_sd,:),ST_Pe(idx_sd,:)); oT=linO(TV_Ps(idx_sd,:),TV_Pe(idx_sd,:));
    sT=10*log10(ST_Ps(idx_sd,:)./ST_Pe(idx_sd,:));
    tT=10*log10(TV_Ps(idx_sd,:)./TV_Pe(idx_sd,:));
    fprintf('【音源%d】\n',idx_sd);
    fprintf('  静的(snapshot1)   : 線形平均 %6.2f / 中央値 %6.2f dB   失敗率 %3.0f%%   BER %.4f\n', ...
            oS, median(sT,'omitnan'), mean(sT<10|isnan(sT))*100, mB(ST_B(idx_sd,:)));
    fprintf('  時変(ドップラー) : 線形平均 %6.2f / 中央値 %6.2f dB   失敗率 %3.0f%%   BER %.4f   （劣化 %+.2f dB）\n\n', ...
            oT, median(tT,'omitnan'), mean(tT<10|isnan(tT))*100, mB(TV_B(idx_sd,:)), oT-oS);
end

%% ==== 星座図 ====
qpsk=[1+1i,1-1i,-1+1i,-1-1i];
figure('Name','時変チャネル版: 静的 vs 時変','Position',[100 100 300*Nsd 560]);
for idx_sd=1:Nsd
    oS=linO(ST_Ps(idx_sd,:),ST_Pe(idx_sd,:)); oT=linO(TV_Ps(idx_sd,:),TV_Pe(idx_sd,:));
    vS=ST_syms{idx_sd}; if isempty(vS)||~any(isfinite(vS)), vS=zeros(4,1); end
    vT=TV_syms{idx_sd}; if isempty(vT)||~any(isfinite(vT)), vT=zeros(4,1); end
    subplot(2,Nsd,idx_sd);
    plot(real(vS),imag(vS),'.','Color',[0.1 0.6 0.9],'MarkerSize',3); hold on;
    plot(real(qpsk),imag(qpsk),'r+','MarkerSize',12,'LineWidth',2); axis equal; grid on;
    lim=safelim(vS); xlim([-lim lim]); ylim([-lim lim]);
    title(sprintf('静的 音源%d\nOSNR=%.1fdB',idx_sd,oS));
    subplot(2,Nsd,Nsd+idx_sd);
    plot(real(vT),imag(vT),'.','Color',[0.85 0.3 0.1],'MarkerSize',3); hold on;
    plot(real(qpsk),imag(qpsk),'r+','MarkerSize',12,'LineWidth',2); axis equal; grid on;
    lim=safelim(vT); xlim([-lim lim]); ylim([-lim lim]);
    title(sprintf('時変(ドップラー) 音源%d\nOSNR=%.1fdB',idx_sd,oT));
end
sgtitle(sprintf('時変チャネル版ver1（f_{wave}=%.2fHz, K=%d, Nf=%d/ch, Nb=%d）',f_wave,K_SUB,Nf,Nb),'FontWeight','bold');

%% ================= ローカル関数 =================
function rx = tv_channel(tx, Hc, Fs, T_wave, Lblk, n_snap)
% 時変畳み込み（ブロック定常近似＋オーバーラップ加算）
%   各ブロック中心時刻の波位相で CIR を補間し、そのCIRでブロックを畳む。
    Ltx=length(tx); Lh=size(Hc,1); rx=zeros(Ltx+Lh-1,1);
    nblk=ceil(Ltx/Lblk);
    for m=1:nblk
        a=(m-1)*Lblk+1; b=min(m*Lblk,Ltx);
        tc=((a+b)/2)/Fs;                         % ブロック中心時刻 [s]
        phi=2*pi*mod(tc,T_wave)/T_wave;          % その時刻の波位相
        h=cir_at_phase(Hc,phi,n_snap);           % CIR補間
        c=fftconv(tx(a:b),h);
        rx(a:a+length(c)-1)=rx(a:a+length(c)-1)+c;
    end
end

function h = cir_at_phase(Hc, phi, n_snap)
% スナップショット列 Hc(:,1..n_snap) を位相 phi で周期線形補間
    f=mod(phi,2*pi)/(2*pi)*n_snap;               % [0, n_snap)
    i0=floor(f); frac=f-i0;
    a=mod(i0,n_snap)+1; b=mod(i0+1,n_snap)+1;
    h=(1-frac)*Hc(:,a)+frac*Hc(:,b);
end

function [comp,psig,perr,ber]=process_beams(sub_ptr,sub_ac,ch_sub,K_SUB,Sps,numSymbols,symbols, ...
                                            gd,Nf,Nb,lambda,delta,Ntrain,mu_pll,freeze,ti)
% サブアレイPTR出力からビーム抽出 → MC-DFE → 評価（ver5と同じ処理）
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
    X=X(:,valid); Kv=size(X,2);
    if Kv==0, comp=zeros(numSymbols,1); psig=1; perr=1e3; ber=0.5; return; end
    out=mc_dfe(X,symbols,Nf,Nb,lambda,delta,Ntrain,mu_pll,freeze);
    if all(isfinite(out))
        [~,ber,comp,psig,perr]=eval_region(symbols,out,ti);
    else
        comp=zeros(numSymbols,1); psig=1; perr=1; ber=0.5;
    end
end

function L = safelim(v)
    a=abs([real(v(:));imag(v(:))]); a=a(isfinite(a));
    if isempty(a), L=3; else, L=max(3,ceil(max(a))); end
end

function h = build_cir_frac(amp, cs, Lh, L)
    h=zeros(Lh,1);
    for p=1:numel(amp)
        c=cs(p); m=(floor(c)-L:floor(c)+L).';
        ok=m>=0 & m<=Lh-1; m=m(ok); if isempty(m), continue; end
        x=m-c; w=0.5*(1+cos(pi*x/(L+1)));
        h(m+1)=h(m+1)+amp(p)*(sinc(x).*w);
    end
end

function y = fftconv(a, b)
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

function out = mc_dfe(X, symbols, Nf, Nb, lambda, delta, Ntrain, mu_pll, freeze)
% マルチチャネルDFE（RLS, ver5=学習後凍結）
    if nargin<9, freeze=true; end
    [N,K]=size(X); Ntap=K*Nf+Nb; w=zeros(Ntap,1);
    for k=1:K, w((k-1)*Nf+1)=1/K; end
    P=(1/delta)*eye(Ntap); ff=zeros(Nf,K); fb=zeros(Nb,1); theta=0; out=zeros(N,1);
    for n=1:N
        ff=[X(n,:); ff(1:end-1,:)];
        ffvec=reshape(ff,[],1)*exp(-1i*theta);
        u=[ffvec; -fb]; y=w'*u; out(n)=y;
        d=sign(real(y))+1i*sign(imag(y));
        if n<=Ntrain, s_ref=symbols(n); fb_new=symbols(n); else, s_ref=d; fb_new=d; end
        if (n<=Ntrain) || ~freeze
            e=s_ref-y; Pu=P*u; kap=Pu/(lambda+u'*Pu);
            w=w+kap*conj(e); P=(P-kap*(Pu'))/lambda;
        end
        phi=imag(y*conj(s_ref)); theta=theta+mu_pll*phi;
        fb=[fb_new; fb(1:end-1)];
    end
end

function [osnr, ber, comp, psig, perr] = eval_region(symbols, raw, idx)
    s=symbols(idx); r=raw(idx); a=(s'*r)/(s'*s); c=r/a;
    psig=mean(abs(s).^2); perr=mean(abs(s-c).^2); osnr=10*log10(psig/perr);
    ber=mean([real(c)>0,imag(c)>0]~=[real(s)>0,imag(s)>0],'all'); comp=raw/a;
end
