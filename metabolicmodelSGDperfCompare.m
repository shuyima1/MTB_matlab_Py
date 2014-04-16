load('Beste7H9modJ.mat')
load('BesteModels.mat', 'BesteVitro')
load('Jamshidimodels.mat', 'jamshidiBiGG')
load('BesteMediaModels.mat', 'BesteExchange')

[iaex iaex ibex] = intersect(BesteExchange(:,1),BesteVitro.rxns);
[iaex2 iaex2 ibex2] = intersect(BesteExchange(:,1),Beste7H9Jaa2.rxns);

% model modification
Bxt = [BesteExchange(iaex2,1) num2cell([BesteVitro.ub(ibex) Beste7H9Jaa2.ub(ibex2)])];

% for 7H9 media
BV2 = BesteVitro;
BV2.ub(ibex) = Beste7H9Jaa2.ub(ibex2);
BV2.c(724) = 1;

% for Griffin media
BVG = BesteVitro;
BVG.ub(ibex) = BesteGriffinJ.ub(ibex2);
BVG.c(724) = 1;

Bxt = [BesteExchange(iaex2,1:2) num2cell([BesteVitro.ub(ibex) Beste7H9Jaa2.ub(ibex2) BV2.ub(ibex) BesteGriffinJ.ub(ibex2) BVG.ub(ibex2)])];

% Jamshidi model modifications
JamshidiExchange = strfind(jamshidiBiGG.rxns,'EX_');
JamshidiExchange = cellfun(@(x) ~isempty(x),JamshidiExchange);
Jxt = [jamshidiBiGG.rxns(JamshidiExchange) printRxnFormula(jamshidiBiGG,Jxt(:,1)) num2cell([jamshidiBiGG.ub(JamshidiExchange) jamshidiBiGG.lb(JamshidiExchange)])];

% for 7H9 media
J7 = changeRxnBounds(jamshidiBiGG,{'EX_leu_L(e)';'EX_lys_L(e)';'EX_met_L(e)';'EX_ile_L(e)';'EX_his_L(e)';'EX_gly(e)';'EX_arg_L(e)';'EX_asn_L(e)';'EX_asp_L(e)';'EX_ala_L(e)';'EX_pro_L(e)';'EX_ser_L(e)';'EX_thr_L(e)';'EX_trp_L(e)';'EX_tyr_L(e)';'EX_val_L(e)'},-10^-4,'l');
J7 = changeRxnBounds(J7,'EX_mg2(e)',-1,'l');
J7 = changeRxnBounds(J7,'EX_h2co3(e)',-1,'l');

% for Griffin media
JG = changeRxnBounds(jamshidiBiGG,'EX_h2co3(e)',-1,'l');
JG = changeRxnBounds(JG,'EX_na1(e)',0,'l');
JG = changeRxnBounds(JG,'EX_glu_L(e)',0,'l');
JG = changeRxnBounds(JG,'EX_glc(e)',0,'l');
JG = changeRxnBounds(JG,'EX_cu2(e)',0,'l');
JG = changeRxnBounds(JG,'EX_cl(e)',0,'l');
JG = changeRxnBounds(JG,'EX_ca2(e)',0,'l');
JG = changeRxnBounds(JG,'EX_mg2(e)',-1,'l');
Jxt = [jamshidiBiGG.rxns(JamshidiExchange) printRxnFormula(jamshidiBiGG,Jxt(:,1)) num2cell([jamshidiBiGG.lb(JamshidiExchange) J7.lb(JamshidiExchange)])];

[grRatioJamGriff grKOJamGriff grWTJamGriff] = singleGeneDeletion(JG);
[grRatioJam7H9 grKOJam7H9 grWTJam7H9] = singleGeneDeletion(J7);
[grRatioBV7H9 grKOBV7H9 grWTBV7H9] = singleGeneDeletion(BV2);
[grRatioBVG grKOBVG grWTBVG] = singleGeneDeletion(BVG);
[grRatioBJG grKOBJG grWTBJG] = singleGeneDeletion(BesteGriffinJ);
[grRatioBJ7H9 grKOBJ7H9 grWTBJ7H9] = singleGeneDeletion(Beste7H9Jaa2);

[cJ iaJ ibJ] = intersect(J7.genes,GriffinData(:,1),'stable');
[cBV iaBV ibBV] = intersect(BV2.genes,GriffinData(:,1),'stable');
[cBJ iaBJ ibBJ] = intersect(Beste7H9Jaa2.genes,GriffinData(:,1),'stable');
cJ(:,2:3) = GriffinData(ibJ,7:8);
cJ(:,4) = num2cell(grRatioJam7H9(iaJ));
cJ(:,5) = num2cell(grRatioJamGriff(iaJ));
cBV(:,2:3) = GriffinData(ibBV,7:8);
cBV(:,4) = num2cell(grRatioBV7H9(iaBV));
cBV(:,5) = num2cell(grRatioBVG(iaBV));
cBJ(:,2:3) = GriffinData(ibBJ,7:8);
cBJ(:,4) = num2cell(grRatioBJ7H9(iaBJ));
cBJ(:,5) = num2cell(grRatioBJG(iaBJ));

threshrange = [0:0.05:1]';
perfJam7H9 = zeros(numel(threshrange),2);
perfBV7H9 = zeros(numel(threshrange),2);
perfBJ7H9 = zeros(numel(threshrange),2);
for i = 1:numel(threshrange)
    perfJam7H9(i,1) = sum(cell2mat(cJ(strcmp('essential',cJ(:,3))|strcmp('growth-defect',cJ(:,3)),4)) < threshrange(i))/sum(strcmp('essential',cJ(:,3))|strcmp('growth-defect',cJ(:,3)));
    perfJam7H9(i,2) = sum(cell2mat(cJ(strcmp('non-essential',cJ(:,3)),4)) >= threshrange(i))/sum(strcmp('non-essential',cJ(:,3)));

    perfBV7H9(i,1) = sum(cell2mat(cBV(strcmp('essential',cBV(:,3))|strcmp('growth-defect',cBV(:,3)),4)) < threshrange(i))/sum(strcmp('essential',cBV(:,3))|strcmp('growth-defect',cBV(:,3)));
    perfBV7H9(i,2) = sum(cell2mat(cBV(strcmp('non-essential',cBV(:,3)),4)) >= threshrange(i))/sum(strcmp('non-essential',cBV(:,3)));

    perfBJ7H9(i,1) = sum(cell2mat(cBJ(strcmp('essential',cBJ(:,3))|strcmp('growth-defect',cBJ(:,3)),4)) < threshrange(i))/sum(strcmp('essential',cBJ(:,3))|strcmp('growth-defect',cBJ(:,3)));
    perfBJ7H9(i,2) = sum(cell2mat(cBJ(strcmp('non-essential',cBJ(:,3)),4)) >= threshrange(i))/sum(strcmp('non-essential',cBJ(:,3)));

    perfJam7H9(i,4) = sum(cell2mat(cJ(strcmp('essential',cJ(:,3))|strcmp('growth-defect',cJ(:,3)),4)) < threshrange(i));
    perfJam7H9(i,5) = sum(cell2mat(cJ(strcmp('non-essential',cJ(:,3)),4)) >= threshrange(i));

    perfBV7H9(i,4) = sum(cell2mat(cBV(strcmp('essential',cBV(:,3))|strcmp('growth-defect',cBV(:,3)),4)) < threshrange(i));
    perfBV7H9(i,5) = sum(cell2mat(cBV(strcmp('non-essential',cBV(:,3)),4)) >= threshrange(i));

    perfBJ7H9(i,4) = sum(cell2mat(cBJ(strcmp('essential',cBJ(:,3))|strcmp('growth-defect',cBJ(:,3)),4)) < threshrange(i));
    perfBJ7H9(i,5) = sum(cell2mat(cBJ(strcmp('non-essential',cBJ(:,3)),4)) >= threshrange(i));
end
perfJam7H9(:,3) = perfJam7H9(:,1)/2+perfJam7H9(:,2)/2;
perfBV7H9(:,3) = perfBV7H9(:,1)/2+perfBV7H9(:,2)/2;
perfBJ7H9(:,3) = perfBJ7H9(:,1)/2+perfBJ7H9(:,2)/2;

sensJamGriffin = zeros(numel(threshrange),3);
sensBVGriffin = zeros(numel(threshrange),3);
sensBJGriffin = zeros(numel(threshrange),3);
specJamGriffin = zeros(numel(threshrange),3);
specBVGriffin = zeros(numel(threshrange),3);
specBJGriffin = zeros(numel(threshrange),3);
truerange = [0.05 0.1 0.15]';
for i = 1:numel(truerange)
    for j = 1:numel(threshrange)
        sensJamGriffin(j,i) = sum(cell2mat(cJ(cell2mat(cJ(:,2)) < truerange(i),5)) < threshrange(j))/sum(cell2mat(cJ(:,2)) < truerange(i));
        specJamGriffin(j,i) = sum(cell2mat(cJ(cell2mat(cJ(:,2)) >= truerange(i),5)) >= threshrange(j))/sum(cell2mat(cJ(:,2)) >= truerange(i));
        sensBVGriffin(j,i) = sum(cell2mat(cBV(cell2mat(cBV(:,2)) < truerange(i),5)) < threshrange(j))/sum(cell2mat(cBV(:,2)) < truerange(i));
        specBVGriffin(j,i) = sum(cell2mat(cBV(cell2mat(cBV(:,2)) >= truerange(i),5)) >= threshrange(j))/sum(cell2mat(cBV(:,2)) >= truerange(i));
        sensBJGriffin(j,i) = sum(cell2mat(cBJ(cell2mat(cBJ(:,2)) < truerange(i),5)) < threshrange(j))/sum(cell2mat(cBJ(:,2)) < truerange(i));
        specBJGriffin(j,i) = sum(cell2mat(cBJ(cell2mat(cBJ(:,2)) >= truerange(i),5)) >= threshrange(j))/sum(cell2mat(cBJ(:,2)) >= truerange(i));
    end
end
perfJamGriffin = sensJamGriffin/2+specJamGriffin/2;
perfBVGriffin = sensBVGriffin/2+specBVGriffin/2;
perfBJGriffin = sensBJGriffin/2+specBJGriffin/2;
