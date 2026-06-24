%% ============================================================
%  PTR 診断スクリプト ver3
%  ★サブアレイ（距離グループ）ごとのPSR・OSNRを個別評価
%  目的：通常PTR の OSNR ギャップがクロスグループ・タイミングズレに
%        よるものかを確認し、SCM-PTR の有効性を裏付ける
% ============================================================
[Arr, Pos] = read_arrivals_asc('JamTank_Arr.arr');

Nrr = size(Arr,1);   % 受信距離数
Nrz = size(Arr,2);   % 受信深度数
Nsd = size(Arr,3);   % 音源数
NUM_RX = Nrr * Nrz;

fprintf('ARR構成: 距離%d × 深度%d = %dch, 音源%d\n', Nrr, Nrz, NUM_RX, Nsd);
if isfield(Pos,'r')&&isfield(Pos.r,'r'), fprintf('  距離(m): '); fprintf('%.1f ',Pos.r.r); fprintf('\n'); end
if isfield(Pos,'r')&&isfield(Pos.r,'z'), fprintf('  深度(m): '); fprintf('%.2f ',Pos.r.z); fprintf('\n\n'); end

% チャネル ⇔ (距離,深度) 対応
ch_rr=zeros(NUM_RX,1); ch_rz=zeros(NUM_RX,1); cc=0;
for irr=1:Nrr, for irz=1:Nrz, cc=cc+1; ch_rr(cc)=irr; ch_rz(cc)=irz; end, end

Rs = 500; Sps = 16; Fs = Rs*Sps; numSymbols = 2000;

% 結果格納
PSR_sub  = zeros(Nsd, Nrr);   % サブアレイ(=距離)ごとのPSR
OSNR_sub = zeros(Nsd, Nrr);   % サブアレイごとのOSNR
PSR_all  = zeros(Nsd, 1);     % 全ch合算のPSR
OSNR_all = zeros(Nsd, 1);     % 全ch合算(通常PTR)のOSNR
OSNR_scm = zeros(Nsd, 1);     % SCM-MRC後のOSNR

fprintf('（各音源について、全ch合算PTR / 距離グループ別PTR / SCM-MRC の PSR・OSNR を表示します）\n');

for idx_sd = 1:Nsd

    bits=randi([0 1],2*numSymbols,1);
    symbols=(2*bits(1:2:end)-1)+1i*(2*bits(2:2:end)-1);
    sigPower=mean(abs(symbols).^2);
    txSignal=zeros(numSymbols*Sps,1);
    for i=1:numSymbols, txSignal((i-1)*Sps+1:i*Sps)=symbols(i); end

    % 全ch共通の最小遅延
    gmin=inf;
    for ch=1:NUM_RX
        d=real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay); d=d(d>0);
        if ~isempty(d), gmin=min(gmin,min(d)); end
    end

    ptr_buf_len=length(txSignal)*4;
    total_ptr=zeros(ptr_buf_len,1);
    total_ac =zeros(ptr_buf_len,1);
    sub_ptr  =zeros(ptr_buf_len, Nrr);
    sub_ac   =zeros(ptr_buf_len, Nrr);

    for ch=1:NUM_RX
        amp=Arr(ch_rr(ch),ch_rz(ch),idx_sd).A;
        dl =real(Arr(ch_rr(ch),ch_rz(ch),idx_sd).delay);
        v=(dl>0); amp=amp(v); dl=dl(v); dl=dl-gmin;
        v=(dl>=0); amp=amp(v); dl=dl(v);
        if isempty(dl), continue; end

        mds=round(max(dl)*Fs);
        rxBase=zeros(length(txSignal)+mds,1);
        for p=1:length(dl)
            ds=round(dl(p)*Fs);
            rxBase(ds+(1:length(txSignal)))=rxBase(ds+(1:length(txSignal)))+amp(p)*txSignal;
        end
        sp=mean(abs(rxBase).^2); np=sp/10^(20/10);
        rx=rxBase+sqrt(np/2)*(randn(size(rxBase))+1i*randn(size(rxBase)));

        ds2=round(dl*Fs)+1; h=zeros(max(ds2),1); h(ds2)=amp;
        g=conj(flipud(h));
        z=conv(rx,g); hac=conv(h,g);
        Lz=min(length(z),ptr_buf_len); La=min(length(hac),ptr_buf_len);

        total_ptr(1:Lz)=total_ptr(1:Lz)+z(1:Lz);
        total_ac(1:La) =total_ac(1:La)+hac(1:La);

        irr = ch_rr(ch);
        sub_ptr(1:Lz,irr)=sub_ptr(1:Lz,irr)+z(1:Lz);
        sub_ac(1:La,irr) =sub_ac(1:La,irr)+hac(1:La);
    end

    %% 全ch合算 PSR
    [qpk,pk]=max(abs(total_ac));
    lags=(-numSymbols:numSymbols); qsym=zeros(length(lags),1);
    for ii=1:length(lags)
        idx=pk+lags(ii)*Sps;
        if idx>=1&&idx<=length(total_ac), qsym(ii)=total_ac(idx); end
    end
    main=abs(qsym(lags==0))^2; side=sum(abs(qsym(lags~=0)).^2);
    PSR_all(idx_sd)=10*log10(main/side);

    %% 全ch合算 通常PTR OSNR（単一オフセット）
    off=pk+round(Sps/2)-1;
    ps=total_ptr(off:Sps:off+Sps*(numSymbols-1));
    a=(symbols'*ps)/(symbols'*symbols); pc=ps/a;
    OSNR_all(idx_sd)=10*log10(sigPower/mean(abs(symbols-pc).^2));

    %% サブアレイ（距離グループ）ごとの個別 PSR・OSNR
    sub_syms = cell(Nrr,1);
    sub_err  = zeros(Nrr,1);
    for irr=1:Nrr
        [qpk_s,pk_s]=max(abs(sub_ac(:,irr)));
        if qpk_s==0, PSR_sub(idx_sd,irr)=NaN; OSNR_sub(idx_sd,irr)=NaN; continue; end

        % PSR
        qsym_s=zeros(length(lags),1);
        for ii=1:length(lags)
            idx2=pk_s+lags(ii)*Sps;
            if idx2>=1&&idx2<=ptr_buf_len, qsym_s(ii)=sub_ac(idx2,irr); end
        end
        m_s=abs(qsym_s(lags==0))^2; sl_s=sum(abs(qsym_s(lags~=0)).^2);
        PSR_sub(idx_sd,irr)=10*log10(m_s/sl_s);

        % OSNR（正しいオフセット）
        off_s=pk_s+round(Sps/2)-1;
        if off_s+Sps*(numSymbols-1) > ptr_buf_len
            OSNR_sub(idx_sd,irr)=NaN; continue;
        end
        ps_s=sub_ptr(off_s:Sps:off_s+Sps*(numSymbols-1), irr);
        a_s=(symbols'*ps_s)/(symbols'*symbols); pc_s=ps_s/a_s;
        err_s=mean(abs(symbols-pc_s).^2);
        OSNR_sub(idx_sd,irr)=10*log10(sigPower/err_s);
        sub_syms{irr}=pc_s; sub_err(irr)=err_s;
    end

    %% SCM-MRC合成
    scm=zeros(numSymbols,1); wsum=0;
    for irr=1:Nrr
        if isempty(sub_syms{irr})||isnan(sub_err(irr))||sub_err(irr)==0, continue; end
        w=1/sub_err(irr);
        scm=scm+w*sub_syms{irr}; wsum=wsum+w;
    end
    if wsum>0
        scm=scm/wsum;
        a_scm=(symbols'*scm)/(symbols'*symbols); scmC=scm/a_scm;
        OSNR_scm(idx_sd)=10*log10(sigPower/mean(abs(symbols-scmC).^2));
    else
        OSNR_scm(idx_sd)=NaN;
    end

    %% 結果表示
    fprintf('\n=== 音源%d ===\n', idx_sd);
    fprintf('  全ch合算  PSR     : %6.2f dB\n', PSR_all(idx_sd));
    fprintf('  全ch合算  OSNR    : %6.2f dB  ← 通常PTR（タイミングズレで劣化）\n', OSNR_all(idx_sd));
    for irr=1:Nrr
        if isfield(Pos,'r')&&isfield(Pos.r,'r')
            dist_label=sprintf('%.0fm', Pos.r.r(irr));
        else
            dist_label=sprintf('距離%d', irr);
        end
        fprintf('  %s群(%dch) PSR  : %6.2f dB\n', dist_label, Nrz, PSR_sub(idx_sd,irr));
        fprintf('  %s群(%dch) OSNR : %6.2f dB\n', dist_label, Nrz, OSNR_sub(idx_sd,irr));
    end
    fprintf('  SCM-MRC後 OSNR  : %6.2f dB  ← 各群を正確なタイミングで合成\n', OSNR_scm(idx_sd));
    fprintf('  ギャップ(SCM vs 全ch PSR): %6.2f dB\n', PSR_all(idx_sd)-OSNR_scm(idx_sd));
end

%% ===== 可視化 =====
figure('Name','PSR vs OSNR 詳細比較','Position',[50 50 220*(Nrr+2) 450]);

% PSR比較（全ch vs 各サブアレイ）
ax1=subplot(1,3,1);
data_psr=[PSR_all, PSR_sub];
bar(data_psr);
sub_names = arrayfun(@(k)sprintf('距離%d群',k),1:Nrr,'UniformOutput',false);
legend_labels=[{'全ch(合算)'}, sub_names];
legend(legend_labels,'Location','best');
set(gca,'XTickLabel',arrayfun(@(k)sprintf('音源%d',k),1:Nsd,'UniformOutput',false));
ylabel('PSR (dB)'); grid on; title('q関数 PSR比較');

% OSNR比較（通常PTR vs 各サブアレイ vs SCM）
ax2=subplot(1,3,2);
data_osnr=[OSNR_all, OSNR_sub, OSNR_scm];
bar(data_osnr);
legend_labels2=[{'通常PTR(単一offset)'}, sub_names, {'SCM-MRC'}];
legend(legend_labels2,'Location','best');
set(gca,'XTickLabel',arrayfun(@(k)sprintf('音源%d',k),1:Nsd,'UniformOutput',false));
ylabel('OSNR (dB)'); grid on; title('OSNR比較（通常PTR / サブアレイ別 / SCM）');

% PSR-OSNR ギャップ
ax3=subplot(1,3,3);
gap_ptr = PSR_all - OSNR_all;
gap_scm = PSR_all - OSNR_scm;
bar([gap_ptr, gap_scm]);
legend({'通常PTR ギャップ','SCM ギャップ'},'Location','best');
set(gca,'XTickLabel',arrayfun(@(k)sprintf('音源%d',k),1:Nsd,'UniformOutput',false));
ylabel('PSR - OSNR (dB)'); grid on; title('ギャップ（小いほど効率的）');
yline(0,'r--');

sgtitle('PTR診断 ver3：サブアレイ別分解', 'FontSize', 12, 'FontWeight', 'bold');
