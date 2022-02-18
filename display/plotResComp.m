%% plot compare
close all
clear
set(groot, 'DefaultAxesFontSize', 16)
set(groot, 'DefaultLineLineWidth', 2)
plotWidth = 1200;
plotHeight = 800;


%Select folder
selpath = uigetdir('res/');
nameFiles = dir([selpath '/' '*_output.mat']);

markerList = ['+','x','o','*','s','d','h'];
lineList = {'-','--',':','-.'};

col = [ 0.901960784313726	0.098039215686275	0.294117647058823;
        0.235294117647059	0.705882352941176	0.294117647058823;
        1                   0.882352941176471	0.098039215686275;
        0                   0.509803921568627	0.784313725490196;
        0.96078431372549    0.509803921568627	0.188235294117647;
        0.568627450980392   0.117647058823529	0.705882352941176;
        0.274509803921569   0.941176470588235	0.941176470588235;
        0.941176470588235   0.196078431372549	0.901960784313726;
        0.823529411764706   0.96078431372549	0.235294117647059;
        0.980392156862745   0.745098039215686	0.831372549019608;
        0                   0.501960784313726	0.501960784313726;
        0.862745098039216   0.745098039215686	1;
        0.666666666666667   0.431372549019608   0.156862745098039;
        1                   0.980392156862745	0.784313725490196;
        0.501960784313726   0                   0;
        0.666666666666667   1                   0.764705882352941;
        0.501960784313726   0.501960784313726	0;
        1                   0.843137254901961	0.705882352941176;
        0                   0                   0.501960784313726;
        0.501960784313726   0.501960784313726	0.501960784313726;
        1                   1                   1]; % from https://sashamaps.net/docs/resources/20-colors/

        


%% iterate on various sol 

nbIt = length(nameFiles);
percMIPS = 40:20:160;

for it = 1:nbIt
    infoFile = nameFiles(it);
    load([infoFile.folder '/' infoFile.name]);

    %init variable for all files 
    if it==1
        % Some parameteres
        t0   =  PAR.sim.startTime.value;  % initial time
        tf   =  PAR.sim.endTime.value; % end time
        N    =  PAR.sim.step.value; % number of dsicretization steps on the horizon length
        dt   =  (tf-t0)/N;
        k    =  1:N;
        k_no1=  2:N;
        t    =  t0:dt:tf;

        delta = 0.001;
        E     = 10000;

        nbControlItem   = PAR.arch.algos.nb+PAR.arch.sensors.nb+1;
        sizeVar         = N*[PAR.arch.algos.nbMode, PAR.arch.sensors.nbMode, PAR.arch.OBC.nbMode, ... % primary variable u
                            PAR.arch.typeDataOut.nb, PAR.arch.algos.nb, PAR.arch.algos.nb, PAR.arch.algos.nb, PAR.arch.memories.nb, ... % secondary variable x
                            PAR.arch.typeDataOut.nb,  PAR.arch.algos.nb,  PAR.arch.algos.nb];
        nbVariable      = sum(sizeVar);
        nbCtrlVariable  = sum(sizeVar(1:nbControlItem));
        
        nbTypeDataOut   = PAR.arch.typeDataOut.nb;
        nbMemories      = PAR.arch.memories.nb;
        maxPower        = PAR.powerGen.offset.value+PAR.powerGen.amplitude.value;
        maxProcRes      = max(PAR.arch.OBC.procRess.value);
        maxDataRate     = max(PAR.arch.memories.maxDataRate.value);
        
        accData         = zeros(nbTypeDataOut*nbIt,N);
        legendTxtAcc    = cell(nbTypeDataOut*nbIt,1);
        
        memUsage        = zeros(nbMemories*nbIt,N);
        legendTxtMem    = cell(nbMemories*nbIt,1);
        
        powerLeft       = zeros(nbIt,N);
        legendTxtPower  = cell(nbIt,1);
        
        procRessLeft    = zeros(nbIt,N);
        legendTxtProcRes= cell(nbIt,1);
        
        memeDataRate    = zeros(nbMemories*nbIt,N);
        legendTxtDR     = cell(nbMemories*nbIt,1);
        
    end
    
    %scenarioName = ['Amp ' num2str(PAR.powerGen.amplitude.value) 'W'   ];
    %   ['Amp ' num2str(PAR.powerGen.amplitude.value) 'W'   ];
    % scenarioName =   ['Off ' num2str(PAR.powerGen.offset.value) 'W'   ];
     scenarioName =     ['Period ' num2str(PAR.powerGen.period.value) 'Steps'   ];
    %   ['Alpha ' num2str(PAR.arch.typeDataOut.alpha(1)*1000) '*10⁻³'   ];
    % scenarioName =   ['MIPS ' num2str(percMIPS(it)) '%'   ];
    %   ['Weight ' num2str(PAR.objFuncWeight.value{1}(1)) ' & ' num2str(PAR.objFuncWeight.value{1}(2))  ];
    
    if maxProcRes < max(PAR.arch.OBC.procRess.value)
        maxProcRes = max(PAR.arch.OBC.procRess.value);
    end
    
    if maxPower < PAR.powerGen.offset.value+PAR.powerGen.amplitude.value
        maxPower        = PAR.powerGen.offset.value+PAR.powerGen.amplitude.value;
    end
    
    % get result per file 
    
    %Acc Data
    idData = nbCtrlVariable+1;
    for i =1:nbTypeDataOut
        legendTxtAcc{(it-1)*PAR.arch.typeDataOut.nb+i} = ['Data Type ' num2str(i) ' - ' scenarioName ];
        accData((it-1)*PAR.arch.typeDataOut.nb+i,:) = reshape(result.x(idData:idData+N-1),[1,N]);
        idData = idData + N;
    end
    
    % Memory storage
    idData = sum(sizeVar(1:nbControlItem+4))+1;
    for i =1:nbMemories
        legendTxtMem{(it-1)*nbMemories+i} = ['Data Type ' num2str(i) ' - ' scenarioName];
        memUsage((it-1)*nbMemories+i,:) = reshape(result.x(idData:idData+N-1),[1,N]);
        idData = idData + N;
    end
    
    %Power
    powerLeft(it,:) = PAR.powerGen.offset.value+PAR.powerGen.amplitude.value*cos(2*pi*k/PAR.powerGen.period.value);
    idData =sum(N*PAR.arch.algos.nbMode) + 1;
    for i =1:PAR.arch.sensors.nb
        nbMode = PAR.arch.sensors.nbMode(i);
        powerLeft(it,:) = powerLeft(it,:) - (reshape(result.x(idData:idData +N*nbMode-1),[nbMode,N])'*PAR.arch.sensors.power.value{i}')';
        idData = idData + N*nbMode;
    end
    
    nbMode = PAR.arch.OBC.nbMode;
    powerLeft(it,:) = powerLeft(it,:) - (reshape(result.x(idData:idData +N*nbMode-1),[nbMode,N])'*PAR.arch.OBC.power.value)';
    legendTxtPower{it} = ['Power Left - ' scenarioName];
    
    % Process ress
    idData =1;
    for i =1:PAR.arch.algos.nb
    
        nbMode = PAR.arch.algos.nbMode(i);
        procRessLeft(it,:) = procRessLeft(it,:) - (reshape(result.x(idData:idData+N*nbMode-1),[nbMode,N])'*PAR.arch.algos.procRess.value{i}')';
        idData = idData + N*nbMode;
    end
    idData =sum(N*[PAR.arch.algos.nbMode  PAR.arch.sensors.nbMode]) + 1;
    nbMode = PAR.arch.OBC.nbMode;
    procRessLeft(it,:) = procRessLeft(it,:) + (reshape(result.x(idData:idData +N*nbMode-1),[nbMode,N])'*PAR.arch.OBC.procRess.value)';
    legendTxtProcRes{it} = ['Proc. Res. Left - ' scenarioName];
    

    % Data Rate
    idData =sum(N*PAR.arch.algos.nbMode) + 1;
    
    for i =1:PAR.arch.sensors.nb
        nbMode = PAR.arch.sensors.nbMode(i);
        memLink = PAR.arch.sensors.memoryLink(i);
        
        memeDataRate((it-1)*nbMemories+memLink,:) = memeDataRate((it-1)*nbMemories+memLink,:)+ (reshape(result.x(idData:idData +N*nbMode-1),[nbMode,N])'*PAR.arch.sensors.dataRate.value{i}')';
        
        idData = idData + N*nbMode;
    end
    
    for i =1:nbMemories
        legendTxtDR{(it-1)*nbMemories+i} = ['Data Rate @ Mem ' num2str(i) ' - ' scenarioName];
    end
    
    clear PAR result;
end

%% Plot

markerListStep = 5:5:N;
markerListStepShift = 0:5:N;
%Acc Data
figure(1)
hold on
for it = 1:nbIt
    for i =1:nbTypeDataOut
        plot(accData((it-1)*nbTypeDataOut+i,:),[lineList{i} markerList(i)], 'MarkerIndices', markerListStepShift+it,'Color',col(it,:))
    end
end
title("Accuracy")
set(gcf, 'Position',  [100, 100, plotWidth, plotHeight])
ylabel('Accuracy [0;1]');
xlabel('Step');
ylim([0;1]);
grid on
legend(legendTxtAcc,'Location','best')
hold off



% Memory storage
figure(2)
hold on
for it = 1:nbIt
    for i =1:nbMemories
        plot(memUsage((it-1)*nbMemories+i,:),[lineList{i} markerList(i)], 'MarkerIndices', markerListStepShift+it,'Color',col(it,:))
    end
end
title("Memory usage")
set(gcf, 'Position',  [100, 100, plotWidth, plotHeight])
ylabel('Memeory usage [Mb]');
xlabel('Step');
%ylim([0;max(PAR.arch.memories.maxStorage.value)]);
grid on
legend(legendTxtMem,'Location','best')
hold off

% Power
figure(3)
hold on
newSeq = 1:nbIt ;%  [1,4,7];%
for it =  newSeq
    %plot(powerLeft(it,:),['-' markerList(it)], 'MarkerIndices', markerListStepShift+it,'Color',col(it,:))
    plot(powerLeft(it,:),lineList{mod(it-1,4)+1},'Color',col(it,:))
end
title("Power")
set(gcf, 'Position',  [100, 100, plotWidth, plotHeight])
ylabel('Power [W]');
xlabel('Step');
ylim([0;maxPower*0.8]);
grid on
legend(legendTxtPower(newSeq),'Location','best')
hold off

% Processing Res
figure(4)
hold on
for it = 1:nbIt
    plot(procRessLeft(it,:),[lineList{mod(it-1,4)+1} markerList(it)], 'MarkerIndices',  markerListStepShift+it,'Color',col(it,:))
end
title("Processing resources")
set(gcf, 'Position',  [100, 100, plotWidth, plotHeight])
ylabel('Processing resources [MIPS]');
xlabel('Step');
ylim([0;maxProcRes*1.1]);
grid on
legend(legendTxtProcRes,'Location','best')
hold off

% Data rate @ mem
figure(5)

hold on
for it = 1:nbIt
    for i =1:nbMemories
        plot(memeDataRate((it-1)*nbMemories+i,:),['-' markerList(it)], 'MarkerIndices',  markerListStepShift+it,'Color',col(it,:))
    end
end
title("Data rate at Memory")
set(gcf, 'Position',  [100, 100, plotWidth, plotHeight])
ylabel('Data rate [Mb/s]');
xlabel('Step');
ylim([0;maxDataRate*1.1]);
grid on
legend(legendTxtDR,'Location','best')
hold off


    




