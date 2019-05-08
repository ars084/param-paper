clc
clear all
close all

addpath('/Volumes/McVeighLab/projects/Anderson/PV_loops')
paramspath = '/Users/andersonscott/Desktop/Papers/Param Paper/';
cd(paramspath);
Outcomes = readtable('CT and Outcomes.xlsx');
CTparam = readtable('LVAD RV CT Study Overview.xlsx','Sheet',1);
CathParam = readtable('LVAD RV CT Study Overview.xlsx','Sheet',3);
EchoParam = readtable('LVAD RV CT Study Overview.xlsx','Sheet',2);
CTparam.RVESV(12) = 248;
CTparam.RVSV{12} = '92';
EchoParam.RVFunction{23} = 'normal';

%%
% collect patients who were included
n = 1;
m = 1;
MACE = table;
noMACE = table;
for i = 1:size(Outcomes,1)
    if strcmp(Outcomes.EXCLUDED{i},'no')
        if strcmp(Outcomes.MACE{i},'true')
            MACE.patient(n) = i;
            MACE.PAPI(n) = CathParam.PAPI(i);
            MACE.RVSWI(n) = CathParam.RVSWI(i);
            MACE.Coupling(n) = Outcomes.Volume_Coupling(i);
            MACE.SizeRV(n) = EchoParam.LVSize_numeric_(i);
            MACE.FunctionRV(n) = EchoParam.LVFunc_numeric_(i);
            MACE.BSA(n) = CTparam.BSA_m2_Mosteller(i);
            n = n +1;
        elseif strcmp(Outcomes.MACE{i},'false')
            noMACE.patient(m) = i;
            noMACE.PAPI(m) = CathParam.PAPI(i);
            noMACE.RVSWI(m) = CathParam.RVSWI(i);
            noMACE.Coupling(m) = Outcomes.Volume_Coupling(i);
            noMACE.SizeRV(m) = EchoParam.LVSize_numeric_(i);
            noMACE.FunctionRV(m) = EchoParam.LVFunc_numeric_(i);
            noMACE.BSA(m) = CTparam.BSA_m2_Mosteller(i);
            m = m +1;
        end
    end
end
%% Get the RVSWI that's CT and Cath Derived

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


