% Anderson Scott
% 4/28/2019
% Code to make figures comparing imaging modalities

clc
clear all
close all

addpath('/Volumes/McVeighLab/projects/Anderson/PV_loops')
paramspath = '/Users/andersonscott/Desktop/Papers/Param Paper/';
cd(paramspath);
CTparam = readtable('LVAD RV CT Study Overview.xlsx','Sheet',1);
CathParam = readtable('LVAD RV CT Study Overview.xlsx','Sheet',3);
EchoParam = readtable('LVAD RV CT Study Overview.xlsx','Sheet',2);
CTparam.RVESV(12) = 248;
CTparam.RVSV{12} = '92';
EchoParam.RVFunction{23} = 'normal';
exclusion = { %Not applicable to MACE scoring b/c...
    'CVC1711021537',... %no recent followup
    'CVC1803141050',... %not usable
    'CVC1803261638',... %no recent followup
    'CVC1805221143',... %no recent followup
    'CVC1805221523',... %post vad
    'CVC1807311059',... %not usable
    'CVC1808011047',... %post vad
    'CVC1901171310',... %restricted
    'CVC1904011519',... %<1 month
    'CVC1904011736',... %<1 month
    'CVC1904081533',... %<1 month
    'CVC1904220946'}; %<1 month 
exclusion = cell2table(exclusion);

% This shows how estimating SV = CO/HR gets worse with increased TR
% Set to false if you don't care about that
PlotEffectTR = true;

%% Section 1: RVSWI
cutoff = 7; %days between Cath and CT
pat_indx_ct = find(CathParam.DaysBetweenCathAndCT <= cutoff);

n=1;
for i = 1:length(pat_indx_ct)
    ind_ex = find(strcmp(CTparam.AnonName{pat_indx_ct(i),:},exclusion{:,:}));
    if ind_ex > 0
        list(n) = i;
        n = n+1;
    end
end
pat_indx_ct(list) = [];

% Correlation of CT and Cath-derived SV
CathSV = 1000 * CathParam.CO(pat_indx_ct) ./ CathParam.HR(pat_indx_ct);
CTSV = CTparam.RVSV(pat_indx_ct);
CTSV = str2double(CTSV);
figure(1)
subplot(3,2,1)
plot(CathSV,CTSV,'k.','markersize',20)
text(180,230,'Cath vs CT SV','fontsize',28,'fontweight','bold')
xlabel('Cath SV','fontsize',24)
ylabel('CT SV','fontsize',24)
unity = 200;
hold on 
plot([0 unity],[0 unity])
hold off
subplot(3,2,2)
BlandAltman(CathSV,CTSV)

% Correlation of CT and Cath-derived pulse pressure
EstPP = CathParam.MPAP(pat_indx_ct) - CathParam.RAP(pat_indx_ct);
CathPP = CathParam.RVs(pat_indx_ct) - CathParam.RVd(pat_indx_ct);
subplot(3,2,3)
plot(EstPP, CathPP,'k.','markersize',20)
text(50,70,'MPAP-RAP vs RV PP','fontsize',28,'fontweight','bold')
xlabel('MPAP-RAP','fontsize',24)
ylabel('RV PP','fontsize',24)
unity = 60;
hold on 
plot([0 unity],[0 unity])
hold off
subplot(3,2,4)
BlandAltman(EstPP, CathPP)

% Correlation of CT+RHC and Cath-derived RVSWI
RVSWI_path = '/Volumes/McVeighLab/projects/Anderson/PV_loops/LVAD_Waveforms/';
cd(RVSWI_path);
pats = dir('CVC*');
for i = 1:length(CTparam.AnonName)
    if ~isempty(strfind([pats.name],CTparam.AnonName{i}))
        cd([RVSWI_path,CTparam.AnonName{i}])
        load('RVSW.mat')
        CTparam.RVSW(i) = RVSW;
        cd .. 
    end
end

CTRVSWI = CTparam.RVSW./CTparam.BSA_m2_Mosteller;
subindx = find(CTRVSWI > 0);
CathRVSWI = CathParam.RVSWI(subindx);
CTRVSWI = CTRVSWI(subindx);

subplot(3,2,5)
plot(CathRVSWI, CTRVSWI,'k.','markersize',20)
text(11, 30,'Cath vs Cath-CT RVSWI','fontsize',28,'fontweight','bold')
xlabel('Cath','fontsize',24)
ylabel('Cath-CT','fontsize',24)
unity = 14;
hold on 
plot([0 unity],[0 unity])
hold off
subplot(3,2,6)
BlandAltman(CathRVSWI, CTRVSWI)
%% Section 2: Coupling
% RV end systolic pressure is approximated by MPAP for Ea
mpap = CathParam.MPAP(pat_indx_ct);
rvpp = CathParam.RVs(pat_indx_ct) - CathParam.RVd(pat_indx_ct);
figure(2)
subplot(3,2,1)
plot(mpap, rvpp, 'k.','markersize',20)
text(65, 73,'Pulse Pressure','fontsize',24,'fontweight','bold')
xlabel('MPAP','fontsize',24)
ylabel('RVs - RVd','fontsize',24)
unity = 63;
hold on 
plot([0 unity],[0 unity])
hold off

subplot(3,2,2)
BlandAltman(mpap, rvpp)

subplot(3,2,3)
plot(mpap./CathSV,rvpp./CTSV,'k.','markersize',20)
text(1.8, 2.5,'Ea Estimation','fontsize',24,'fontweight','bold')
xlabel('MPAP / SV_C_a_t_h','fontsize',24)
ylabel('RVs / SV_C_T','fontsize',24)
unity = 2;
hold on 
plot([0 unity],[0 unity])
hold off

subplot(3,2,4)
BlandAltman(mpap./CathSV,rvpp./CTSV)

ESV = CTparam.RVESV(pat_indx_ct);
rvs = CathParam.RVs(pat_indx_ct);
rvd = CathParam.RVd(pat_indx_ct);
VCC = CathSV./ESV;
subplot(3,2,5)
plot(VCC,(CTSV./ESV).*(rvs./(rvs - rvd)),'k.','markersize',20)
text(1.3, 2.4,'Volumetric Coupling','fontsize',24,'fontweight','bold')
xlabel('VCC = SV / ESV','fontsize',24)
ylabel('VCC_C_T','fontsize',24)
unity = 1.5;
hold on 
plot([0 unity],[0 unity])
hold off

subplot(3,2,6)
BlandAltman(VCC,(CTSV./ESV).*(rvs./(rvs - rvd)))

if PlotEffectTR == true
    figure(6)
    TR = EchoParam.TricuspRegur(pat_indx_ct);
    TR = str2double(TR)>=2;
    boxplot(CTSV-CathSV,TR)
    title('Diff in SV as a func of TR', 'fontsize',24)
    xlabel('TR','fontsize',20)
    ylabel('\Delta SV','fontsize',20)
end

%% Section 3: Echocardiography

cutoff = 7; %days between Echo and CT
pat_indx_echo = find(strcmp(EchoParam.W_in7Days_,'yes'));

n=1;
for i = 1:length(pat_indx_echo)
    ind_ex_echo = find(strcmp(CTparam.AnonName{pat_indx_echo(i),:},exclusion{:,:}));
    if ind_ex_echo > 0
        list_echo(n) = i;
        n = n+1;
    end
end
pat_indx_echo(list_echo) = [];

% RV size
RVsize = EchoParam.RVSize_numeric_(pat_indx_echo);
EDV = CTparam.RVEDV(pat_indx_echo);
EDV = str2double(EDV);
BSA = CTparam.BSA_m2_Mosteller(pat_indx_echo);
EDVI = EDV./BSA;
figure(3)
subplot(1,2,1)
boxplot(EDV, RVsize>0)
title('RV Size','fontsize',24)
ylabel('EDV (mL)','fontsize',24)
xticklabels({'Normal','Enlarged'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',20)
% RV function
RVfunc = EchoParam.RVFunction(pat_indx_echo);
RVEF = CTparam.RVEF(pat_indx_echo);
subplot(1,2,2)
boxplot(RVEF,RVfunc)
title('RV Function','fontsize',24)
ylabel('RVEF (%)','fontsize',24)
xticklabels({'Reduced','Normal'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',20)
