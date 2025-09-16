resp = 'Median_Accumbens_Area_Left_SUL_LOG';
C = cal.(resp);
fprintf('alpha=%.4f  beta=%.4f  R2_link=%.3f  n_loo=%d\n', C.alpha, C.beta, C.R2_link, C.n_loo_pairs);
%%
TT = cv.(resp);
z975 = norminv(0.975);
y  = double(TT.y_true);
p  = double(TT.y_pred);
sd = double(TT.sd_pred_link);
ok = isfinite(y) & isfinite(p) & isfinite(sd) & sd>0;

ycal = C.alpha + C.beta.*p(ok);
PIlo = ycal - z975*(C.c*sd(ok));
PIhi = ycal + z975*(C.c*sd(ok));
covPI_link = mean(y(ok) >= PIlo & y(ok) <= PIhi);
fprintf('Coverage PI (link) = %.1f%%\n', 100*covPI_link);

%%

% vytvoř predikční křivku pro mediánové kovariáty (např. Sex=M/F zvlášť)
% můžeš použít uložený PNG, ale tady je rychlý overlay pro M (pokud je Sex v datech):
L = M(resp);
ref = refRowLike(L, Tc, resp);    % měl bys ho už mít z reportu; jinak počítej stejně

ageTrain = double(Tc.Age); ageTrain = ageTrain(isfinite(ageTrain));
age = linspace(prctile(ageTrain,5), prctile(ageTrain,95), 100)';

newT = repmat(ref, numel(age), 1);
newT.Age = age;
% doplň AgeR* nebo cAge* podle toho, co má model (viz report_model.m)
feNames = string(L.Variables.Properties.VariableNames);
useRCS  = any(startsWith(feNames,"AgeR"));
if useRCS
    K = sum(startsWith(feNames,"AgeR")) + 1;
    load('./models/_globals/age_knots.mat','knots');
    B = rcsDesign(age, knots);
    for j=1:size(B,2), newT.(sprintf('AgeR%d',j)) = B(:,j); end
else
    muAge = mean(ageTrain,'omitnan');
    if ismember('cAge', newT.Properties.VariableNames),  newT.cAge  = age - muAge; end
    if ismember('cAge2',newT.Properties.VariableNames), newT.cAge2 = (age - muAge).^2; end
    if ismember('cAge3',newT.Properties.VariableNames), newT.cAge3 = (age - muAge).^3; end
end

[yhat_l, yCI_l] = predict(L, newT, 'Conditional', false, 'Alpha', 0.05);
z975 = norminv(0.975);
SEm  = (yCI_l(:,2)-yCI_l(:,1))/(2*z975);
sd_l = sqrt(max(0, SEm.^2 + L.MSE + addedREvariance(L, newT)));

yhat_cal_l = C.alpha + C.beta.*yhat_l;      % link-kalibrace
tr = @(x) exp(x).*smear;

hold on
plot(age, tr(yhat_cal_l), 'k-', 'LineWidth', 2);
