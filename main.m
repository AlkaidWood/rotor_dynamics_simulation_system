clc
clear
close all

%% Input initial parameters

InitialParameter = inputEssentialParameter();
InitialParameter = inputIntermediateBearing(InitialParameter);

% If you need change some parameters, please change the data in the struct:
% InitialParameter, then use establishModel( ) to get the different models

%% Establish models

% grid{1} = [1 4 1]; manualGrid{2} = [3]; 
grid = 'low';
Parameter = establishModel(InitialParameter,...
                           'gridfineness', grid,...
                           'isPlotModel',  true,...
                           'isPlotMesh',   true);
save('modelParameter','Parameter')                       
%%  Generate the dynamic equation

generateDynamicEquation(Parameter);                      
 
%% Calculate response

TSTART = 0;
TEND = 1;
SAMPLINGFREQUENCY = 105000;
ISPLOTSTATUS = true;
REDUCEINTERVAL = 1;
tic
[q, dq, t, convergenceStr] = calculateResponse(Parameter, [TSTART, TEND], SAMPLINGFREQUENCY,ISPLOTSTATUS,REDUCEINTERVAL);
toc

if ~isempty(convergenceStr)
    fprintf('%s \n', convergenceStr)
end

%save('response','t','q','dq')

%% Post Proccessing

% signalProccessing
tSpan = [0 1];
SwitchFigure.displacement       = true;
SwitchFigure.axisTrajectory     = true;
SwitchFigure.axisTrajectory3d   = false;
SwitchFigure.phase              = false;
SwitchFigure.fftSteady          = true;
SwitchFigure.fftTransient       = false;
SwitchFigure.poincare           = false;
SwitchFigure.saveFig            = true;

signalProcessing(q, dq, t,...
                 Parameter, SwitchFigure, tSpan, SAMPLINGFREQUENCY,...
                 'reduceInterval', REDUCEINTERVAL,...
                 'fftTimeInterval', 0.5,...
                 'fftisPlot3DTransient', false,...
                 'fftSuperpositionRatio', 0.5,...
                 'fftXlim', 400,...
                 'isPlotInA4', true,...
                 'fftSteadyLog', true)


%% Real-time monitor

% plot part of response
close all
% dofNo = 2;
% 
% figure
% plot(t,q(dofNo,:))
% set(gcf,'unit','centimeters','position',[3 12 35 8])
% set(gca,'Position',[.05 .1 .92 .75])
% 
% figure
% plot(t,dq(dofNo,:))
% set(gcf,'unit','centimeters','position',[3 1.5 35 8])
% set(gca,'Position',[.05 .1 .92 .75])
% 
% figure
% plot(q(dofNo*4-3,:), q(dofNo*4-2,:))
% set(gcf,'unit','centimeters','position',[3 3 20 16])
% set(gca,'Position',[.05 .1 .92 .75])
