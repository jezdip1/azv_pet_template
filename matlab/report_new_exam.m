function report_new_exam(Tc, Tnew, resp, M, cal, outdir)
% Tc    : tréninková tabulka (kvůli distribuci)
% Tnew  : table s 1 řádkem (nový pacient)
% resp  : region (např. 'Median_Accumbens_Area_Left_SUL_LOG')
% M,cal : modely a kalibrace
% outdir: adresář pro uložení grafů a JSON

    if ~exist(outdir,'dir'), mkdir(outdir); end
    P = predict_new_exam(Tnew, resp, M, cal);

    % ---- čísla uložit do JSON ----
    fid = fopen(fullfile(outdir,[resp '_prediction.json']),'w');
    fwrite(fid, jsonencode(P,'PrettyPrint',true));
    fclose(fid);

    % ---- distribuce tréninku ----
    isLog = endsWith(resp,'_LOG') || endsWith(resp,'_SUL_LOG');
    smear = cal.(resp).smear;
    if ~isfinite(smear) || smear<=0, smear=1; end
    if isLog
        trainVals = exp(double(Tc.(resp))) * smear;
    else
        trainVals = double(Tc.(resp));
    end

    % ---- nová hodnota (pokud ji máme) ----
    if ismember(resp, Tnew.Properties.VariableNames)
        if isLog
            yobs = exp(double(Tnew.(resp))) * smear;
        else
            yobs = double(Tnew.(resp));
        end
    else
        yobs = NaN;
    end

    % ---- graf: boxplot + pacient ----
    f = figure('Visible','off','Color','w');
    hold on
    boxchart(trainVals,'BoxWidth',0.5,'MarkerStyle','.');
    if isfinite(yobs)
        scatter(1, yobs, 80,'r','filled','MarkerEdgeColor','k');
    end
    % přidej PI
    plot([0.7 1.3],[P.PI_orig(1) P.PI_orig(1)],'b--','LineWidth',1.5);
    plot([0.7 1.3],[P.PI_orig(2) P.PI_orig(2)],'b--','LineWidth',1.5);
    % popisky
    ylabel(resp,'Interpreter','none');
    title(sprintf('%s — Pred=%.3f, z=%.2f, p=%.3f', ...
        strrep(resp,'_','-'), P.Pred_orig, P.z, P.p_two_sided), ...
        'Interpreter','none');
    xlim([0.5 1.5]); set(gca,'XTick',[]);
    saveas(f, fullfile(outdir,[resp '_boxplot.png']));
    close(f);
end
