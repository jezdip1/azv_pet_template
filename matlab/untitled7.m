% vyber pacienta (jeden řádek) – např. podle UNIS:
% one = Tc(find(Tc.UNIS=="SOME_ID",1), :);
one = Tc(1, :);

outdir = fullfile('reports','new_exam_SOME_ID');
S = report_new_patient_all_regions(Tc, one, info, M, cal, outdir);
% S = report_new_patient_all_regions(Ttrain, Tnew, info, M, cal, outdir)

disp(S.summary_table(1:10, {'Region','Observed','Pred','PI_lo','PI_hi','z','isOut95','isOut99'}))
