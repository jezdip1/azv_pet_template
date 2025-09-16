function S = report_new_patient_all_regions(Ttrain, Tnew, info, M, cal, outdir)
% Jedno vyšetření (1 řádek v Tnew) -> report pro všechny regiony v info.responses
% - Ukládá per-region boxploty + PI, JSON a CSV souhrn, a z-score "waterfall" plot.
% - Prezentace v orig SUL (tj. exp(link) * smear).
%
% Vstupy:
%   Ttrain : tréninková tabulka (kvůli distribucím pro boxplot)
%   Tnew   : table s jedním řádkem (nový pacient) se STEJNÝMI názvy kovariátů
%   info   : struct s polem .responses (string/cellstr) – názvy DV (…_SUL_LOG)
%   M      : containers.Map {'resp' -> LinearMixedModel}
%   cal    : struct s poli per region (alpha, beta, c, smear, …)
%   outdir : kam uložit výstupy
%
% Výstup:
%   S: struct se souhrnnou tabulkou a cestami k souborům

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
boxplot_files = strings(nR,1);

% === loop přes regiony ===
for i = 1:nR
    resp = responses(i);
    REG(i) = resp;
    try
        % --- predikce a klasifikace (na link, pak exp+smear -> SUL orig)
        P = predict_new_exam(Tnew, char(resp), M, cal);

        % uložit čísla
        Pred(i) = P.Pred_orig;
        CIlo(i) = P.CI_orig(1); CIhi(i) = P.CI_orig(2);
        PIlo(i) = P.PI_orig(1); PIhi(i) = P.PI_orig(2);
        z(i)    = P.z;          pval(i) = P.p_two_sided;
        is95(i) = P.is_outlier_95; is99(i) = P.is_outlier_99;

        % pozorovaná hodnota (pokud je v Tnew)
        if ismember(resp, Tnew.Properties.VariableNames)
            isLog = endsWith(resp,'_LOG') || endsWith(resp,'_SUL_LOG');
            if isLog
                smear = cal.(char(resp)).smear;
                if ~isfinite(smear) || smear<=0, smear=1; end
                Obs(i) = exp(double(Tnew.(char(resp)))) * smear;
            else
                Obs(i) = double(Tnew.(char(resp)));
            end
        end

        % --- per-region boxplot s pacientem + PI ---
        isLog = endsWith(resp,'_LOG') || endsWith(resp,'_SUL_LOG');
        if isLog
            sm = cal.(char(resp)).smear; if ~isfinite(sm) || sm<=0, sm=1; end
            trainVals = exp(double(Ttrain.(char(resp)))) * sm;
        else
            trainVals = double(Ttrain.(char(resp)));
        end

        % robustní ořez extrémů jen pro kresbu (kvůli ose)
        v = isfinite(trainVals);
        pr = prctile(trainVals(v), [1 99]);
        yl = [max(min(trainVals(v)), pr(1))  min(max(trainVals(v)), pr(2))];
        pad = 0.05*(yl(2)-yl(1)); yl = yl + [-pad pad];

        f = figure('Visible','off','Color','w'); hold on
        % boxchart z vektoru – aby to bylo “jedno políčko”
        boxchart(ones(sum(v),1), trainVals(v), 'BoxWidth',0.4, 'MarkerStyle','.');
        if isfinite(Obs(i))
            scatter(1, Obs(i), 70, 'r', 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.9);
        end
        % predikční interval vodorovnými čarami přes políčko
        plot([0.75 1.25], [PIlo(i) PIlo(i)], 'b--','LineWidth',1.6);
        plot([0.75 1.25], [PIhi(i) PIhi(i)], 'b--','LineWidth',1.6);

        % popisky (bez TeXu)
        ylabel(safe_label(y_label_sul(resp)), 'Interpreter','none');
        title(safe_label(sprintf('%s  |  Pred=%.3g  z=%.2f  p=%.3g', ...
                 strip_suffix(resp), Pred(i), z(i), pval(i))), 'Interpreter','none');
        xlim([0.5 1.5]); set(gca,'XTick',[]);
        ylim(yl); grid on

        fn = fullfile(outdir, [safe_fname(char(resp)) '_boxplot.png']);
        saveas(f, fn); close(f);
        boxplot_files(i) = string(fn);

        % ulož i JSON s kompletním P
        fid = fopen(fullfile(outdir,[safe_fname(char(resp)) '_prediction.json']),'w');
        fwrite(fid, jsonencode(P,'PrettyPrint',true));
        fclose(fid);

    catch ME
        warning('report: %s failed -> %s', char(resp), ME.message);
    end
end

% === souhrnná tabulka (seřazená dle |z| k největším odchylkám) ===
Tsum = table(REG, Obs, Pred, CIlo, CIhi, PIlo, PIhi, z, pval, is95, is99, boxplot_files, ...
    'VariableNames', {'Region','Observed','Pred','CI_lo','CI_hi','PI_lo','PI_hi','z','p','isOut95','isOut99','boxplot_png'});

[~,ord] = sort(abs(Tsum.z), 'descend', 'MissingPlacement','last');
Tsum = Tsum(ord,:);

writetable(Tsum, fullfile(outdir,'new_exam_summary.csv'));
fid = fopen(fullfile(outdir,'new_exam_summary.json'),'w');
fwrite(fid, jsonencode(Tsum, 'PrettyPrint', true)); fclose(fid);

% === waterfall / hřeby: z-scores přes regiony ===
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
    warning('waterfall plot failed: %s', ME.message);
end

% === návratová struktura ===
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
    s = strrep(txt,'_','-');
    s = strrep(s,'^',' ');
end
function s = safe_fname(txt)
    s = regexprep(txt,'[^A-Za-z0-9\-]+','_');
end
function t = safe_tiny(txt)
    t = char(safe_label(txt));
    if numel(t)>22, t = [t(1:20) '…']; end
end
