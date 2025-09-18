function S = report_new_patient_all_regions(Ttrain, Tnew, info, M, cal, outdir)
% Jedno vyšetření (1 řádek v Tnew) -> report pro všechny regiony
% Uloží: per-region boxplot (MODEL) + boxplot (MATCHED), JSON predikci, CSV/JSON sumáře, z-waterfall.

if ~exist(outdir,'dir'), mkdir(outdir); end
responses = string(info.responses(:));
nR = numel(responses);

% Výstupní pole
REG     = strings(nR,1);
Pred    = nan(nR,1); CIlo = nan(nR,1); CIhi = nan(nR,1);
PIlo    = nan(nR,1); PIhi = nan(nR,1);
Obs     = nan(nR,1);
z       = nan(nR,1); pval = nan(nR,1);
is95    = false(nR,1); is99 = false(nR,1);
box_model_files   = strings(nR,1);
box_matched_files = strings(nR,1);

% === loop přes regiony ===
for i = 1:nR
    resp = responses(i);
    REG(i) = resp;
    try
        % --- predikce a klasifikace (link -> orig přes exp*smear) ---
        P = predict_new_exam(Tnew, char(resp), M, cal);

        % uložit čísla
        Pred(i) = P.Pred_orig;
        CIlo(i) = P.CI_orig(1); CIhi(i) = P.CI_orig(2);
        PIlo(i) = P.PI_orig(1); PIhi(i) = P.PI_orig(2);
        z(i)    = P.z;          pval(i) = P.p_two_sided;
        is95(i) = P.is_outlier_95; is99(i) = P.is_outlier_99;

        % pozorovaná hodnota (pokud je v Tnew)
        isLog = endsWith(resp,'_LOG') || endsWith(resp,'_SUL_LOG');
        smear = 1;
        if isLog, smear = cal.(char(resp)).smear; if ~isfinite(smear) || smear<=0, smear=1; end, end
        Obs(i) = NaN;
        if ismember(resp, Tnew.Properties.VariableNames)
            Obs(i) = double(Tnew.(char(resp)));
            if isLog, Obs(i) = exp(Obs(i)) * smear; end
        end

        % ========= A) MODEL-BASED PREDICTIVE BOXPLOT =========
        Q = predictive_box(cal.(char(resp)), P.Pred_link, P.sd_pred_link, isLog);
        % Y-limření (rozumný padding okolo prediktivního 2.5–97.5 %)
        % ylA = [Q.q2p5 Q.q97p5]; padA = 0.08*(ylA(2)-ylA(1)); ylA = ylA + [-padA padA];
        % 
        % f = figure('Visible','off','Color','w'); hold on
        % % „ruční“ box: Q1–Q3, median, whiskery 2.5–97.5 %
        % line([1 1],[Q.q25 Q.q75],'LineWidth',12,'Color',[0.85 0.92 1]);
        % plot([0.85 1.15],[Q.q50 Q.q50],'b-','LineWidth',2);                 % median
        % plot([1 1],[Q.q2p5 Q.q25],'k-','LineWidth',1.2);                    % whisker low
        % plot([1 1],[Q.q75  Q.q97p5],'k-','LineWidth',1.2);                   % whisker high
        % yline(Q.mean_pred_orig, 'b--','LineWidth',1.6);                      % dashed mean pred
        % if isfinite(Obs(i)), plot(1, Obs(i), 'ro','MarkerFaceColor','r'); end
        % xlim([0.5 1.5]); set(gca,'XTick',[]);
        % ylim(ylA); grid on
        valsA = [Q.q2p5, Q.q97p5, Q.q25, Q.q75, Q.q50, Q.mean_pred_orig, ...
         Pred(i), PIlo(i), PIhi(i), Obs(i)];
        ylA = nice_ylim(valsA, 0.08, [0 1]);   % 8% padding, fallback [0 1]
        
        f = figure('Visible','off','Color','w'); hold on
        % „ruční“ box:
        line([1 1],[Q.q25 Q.q75],'LineWidth',12,'Color',[0.85 0.92 1]);
        plot([0.85 1.15],[Q.q50 Q.q50],'b-','LineWidth',2);
        plot([1 1],[Q.q2p5 Q.q25],'k-','LineWidth',1.2);
        plot([1 1],[Q.q75  Q.q97p5],'k-','LineWidth',1.2);
        yline(Q.mean_pred_orig,'b--','LineWidth',1.6);
        if isfinite(Obs(i)), plot(1, Obs(i), 'ro','MarkerFaceColor','r'); end
        xlim([0.5 1.5]); set(gca,'XTick',[]);
        ylim(ylA); grid on
        ylabel(safe_label(y_label_sul(resp)), 'Interpreter','none');
        title(safe_label(sprintf('%s  |  Pred=%.3g  z=%.2f  p=%.3g (MODEL)',...
            strip_suffix(resp), Pred(i), z(i), pval(i))), 'Interpreter','none');
        fnA = fullfile(outdir, [safe_fname(char(resp)) '_box_model.png']);
        saveas(f, fnA); close(f);
        box_model_files(i) = string(fnA);

        % ========= B) MATCHED EMPIRICAL BOXPLOT =========
        idx  = matched_subset(Ttrain, Tnew, 5);    % ±5 let, stejné pohlaví & scanner (pokud jsou)
        vals = nan(sum(idx),1);
        if sum(idx)>0
            if isLog
                tmp  = Ttrain.(char(resp));       % <- krok 1
                vals = exp(double(tmp(idx))) * smear;  % <- krok 2
            else
                tmp  = Ttrain.(char(resp));
                vals = double(tmp(idx));
            end
        end
        f = figure('Visible','off','Color','w'); hold on
        if ~isempty(vals) && any(isfinite(vals))
            boxchart(ones(numel(vals),1), vals, 'BoxWidth',0.4, 'MarkerStyle','.');
        else
            text(1, 0.5, 'No matched data', 'HorizontalAlignment','center');
        end
        yline(Pred(i), 'b--','LineWidth',1.6);
        if isfinite(Obs(i)), scatter(1, Obs(i), 70, 'r','filled','MarkerEdgeColor','k'); end
        xlim([0.5 1.5]); set(gca,'XTick',[]);
        
        % limity z matched dat (pokud jsou), vždy ale přimícháme Pred/PI/Obs
        baseB = [Pred(i), PIlo(i), PIhi(i), Obs(i)];
        if any(isfinite(vals))
            pr = prctile(vals(isfinite(vals)), [1 99]);
            ylB = nice_ylim([pr(:).', baseB], 0.08, [min(vals) max(vals)]);
        else
            ylB = nice_ylim([ylA, baseB], 0.08, ylA);  % fallback na modelové limity
        end
        ylim(ylB); grid on

        ylabel(safe_label(y_label_sul(resp)), 'Interpreter','none');
        title(safe_label(sprintf('%s  |  Pred=%.3g  z=%.2f  p=%.3g (MATCHED)',...
            strip_suffix(resp), Pred(i), z(i), pval(i))), 'Interpreter','none');
        fnB = fullfile(outdir, [safe_fname(char(resp)) '_box_matched.png']);
        saveas(f, fnB); close(f);
        box_matched_files(i) = string(fnB);

        % --- ulož i JSON s kompletním P ---
        fid = fopen(fullfile(outdir,[safe_fname(char(resp)) '_prediction.json']),'w');
        fwrite(fid, jsonencode(P,'PrettyPrint',true)); fclose(fid);

    catch ME
        warning('report: %s failed -> %s', char(resp), ME.message);
    end
end

% === souhrnná tabulka (seřazená dle |z|) ===
Tsum = table(REG, Obs, Pred, CIlo, CIhi, PIlo, PIhi, z, pval, is95, is99, ...
             box_model_files, box_matched_files, ...
    'VariableNames', {'Region','Observed','Pred','CI_lo','CI_hi','PI_lo','PI_hi','z','p','isOut95','isOut99','box_model_png','box_matched_png'});

[~,ord] = sort(abs(Tsum.z), 'descend', 'MissingPlacement','last');
Tsum = Tsum(ord,:);

writetable(Tsum, fullfile(outdir,'new_exam_summary.csv'));
fid = fopen(fullfile(outdir,'new_exam_summary.json'),'w');
fwrite(fid, jsonencode(Tsum, 'PrettyPrint', true)); fclose(fid);

% === z-waterfall ===
try
    f = figure('Visible','off','Color','w'); hold on
    zr = Tsum.z; zr(~isfinite(zr)) = NaN;
    x = 1:height(Tsum);
    stem(x, zr, 'filled','LineWidth',1.1);
    yline(0,'k-'); yline(1.96,'r--'); yline(-1.96,'r--');
    yline(2.576,'r-'); yline(-2.576,'r-');
    xlim([0 height(Tsum)+1]);
    set(gca,'XTick',x,'XTickLabel',arrayfun(@(s)safe_tiny(strip_suffix(Tsum.Region(s))), (1:height(Tsum))', 'uni',0), ...
            'XTickLabelRotation',90);
    ylabel('z-score (link, calibrated)'); grid on
    title('New exam — standardized residuals across regions','Interpreter','none');
    saveas(f, fullfile(outdir,'new_exam_z_waterfall.png')); close(f);
catch ME
    warning('waterfall plot failed: %s', E.message);
end

S = struct('summary_table', Tsum, ...
           'csv', fullfile(outdir,'new_exam_summary.csv'), ...
           'json', fullfile(outdir,'new_exam_summary.json'), ...
           'z_plot', fullfile(outdir,'new_exam_z_waterfall.png'));
end

% ---------- helpers ----------
function s = strip_suffix(resp)
    s = regexprep(char(resp),'_SUL_LOG$','-SUL');
end
function s = y_label_sul(resp)
    s = [regexprep(char(resp),'_SUL_LOG$','-SUL') ' (orig)'];
end
function s = safe_label(txt)
    s = strrep(txt,'_','-'); s = strrep(s,'^',' ');
end
function s = safe_fname(txt)
    s = regexprep(txt,'[^A-Za-z0-9\-]+','_');
end
function t = safe_tiny(txt)
    t = char(safe_label(txt));
    if numel(t)>22, t = [t(1:20) '…']; end
end
function Q = predictive_box(C, yhat_link, sd_pred_link, isLog)
    s = C.c .* sd_pred_link;
    zf = @(p) norminv(p);
    qL = [0.025 0.25 0.5 0.75 0.975];
    q_link = yhat_link + s.*zf(qL);
    if isLog
        smear = C.smear; if ~isfinite(smear) || smear<=0, smear = 1; end
        tr = @(x) exp(x).*smear;
    else, tr = @(x) x; end
    Q.q2p5  = tr(q_link(1)); Q.q25 = tr(q_link(2)); Q.q50 = tr(q_link(3));
    Q.q75   = tr(q_link(4)); Q.q97p5= tr(q_link(5));
    Q.mean_pred_orig = tr(yhat_link);
end
% function idx = matched_subset(T, Tnew, winY)
% % Podmíněná empirická kohorta: ±winY let, stejný sex & scanner (pokud existují)
%     idx = true(height(T),1);
%     if ismember('Age', T.Properties.VariableNames)
%         idx = idx & abs(double(T.Age) - double(Tnew.Age)) <= winY;
%     end
%     if ismember('Sex', T.Properties.VariableNames) && ismember('Sex', Tnew.Properties.VariableNames)
%         idx = idx & (T.Sex == Tnew.Sex);
%     end
%     if ismember('scanner', T.Properties.VariableNames) && ismember('scanner', Tnew.Properties.VariableNames)
%         idx = idx & (T.scanner == Tnew.scanner);
%     end
% end
function idx = matched_subset(Ttrain, Trow, tolYears)
% Vrátí masku řádků Ttrain, které odpovídají Trow:
% - |Age_train - Age_row| <= tolYears (pokud Age existuje a je definováno)
% - stejné pohlaví (pokud je sloupec Sex v obou tabulkách)
% - stejný scanner (pokud existuje v obou — s aliasy)
%
% Když Ttrain je prázdný nebo není table, vrátí log. prázdný vektor.

    % default prázdný výstup
    if ~istable(Ttrain) || height(Ttrain)==0
        idx = false(0,1);
        return;
    end
    if ~istable(Trow) || height(Trow)~=1
        error('matched_subset: Trow must be a 1-row table.');
    end
    if nargin<3 || isempty(tolYears), tolYears = 5; end

    vnT  = string(Ttrain.Properties.VariableNames);
    vnR  = string(Trow.Properties.VariableNames);

    idx = true(height(Ttrain),1);

    % --- Age filtr ---
    if ismember("Age", vnT) && ismember("Age", vnR) && isfinite(double(Trow.Age))
        idx = idx & abs(double(Ttrain.Age) - double(Trow.Age)) <= tolYears;
    end

    % --- Sex filtr (M/F) ---
    if ismember("Sex", vnT) && ismember("Sex", vnR)
        sxR = Trow.Sex;
        if ~iscategorical(sxR), sxR = categorical(string(sxR),["M","F"]); end
        sxT = Ttrain.Sex;
        if ~iscategorical(sxT), sxT = categorical(string(sxT),["M","F"]); end
        idx = idx & (sxT == sxR);
    end

    % --- Scanner filtr (volitelně; použij aliasy) ---
    scanAliases = ["scanner","Scanner","ScannerID","ScannerName","ScannerModel","ScannerType"];
    colT = scanAliases(ismember(scanAliases, vnT));
    colR = scanAliases(ismember(scanAliases, vnR));
    if ~isempty(colT) && ~isempty(colR)
        st = Ttrain.(char(colT(1)));
        sr = Trow.(char(colR(1)));
        if ~iscategorical(st), st = categorical(string(st)); end
        if ~iscategorical(sr), sr = categorical(string(sr)); end
        idx = idx & (st == sr);
    end
end


function yl = nice_ylim(vals, padFrac, fallback)
% Bezpečné Y-lim z vektoru hodnot. Zahrnuje padding, řeší NaN/Inf,
% nulový rozsah i obrácené pořadí.
    v = vals(isfinite(vals));
    if isempty(v)
        yl = fallback;
        if numel(yl)~=2 || ~all(isfinite(yl)) || yl(1)>=yl(2)
            yl = [0 1];  % úplný poslední fallback
        end
        return
    end
    lo = min(v); hi = max(v);
    if lo == hi
        % rozšířit o ±5 % (min. epsilon)
        span = max(abs(lo), 1) * 0.05;
        if span == 0, span = 1e-6; end
        lo = lo - span; hi = hi + span;
    end
    pad = padFrac * (hi - lo);
    yl = [lo - pad, hi + pad];
    if ~all(isfinite(yl)) || yl(1) >= yl(2)
        % pojistka
        yl = [lo hi];
        if yl(1) == yl(2)
            yl = yl + [-1 1]*max(abs(yl(1)),1)*0.05;
        end
    end
end
