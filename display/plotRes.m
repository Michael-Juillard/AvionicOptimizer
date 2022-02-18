%% plot
close all
clear
set(groot, 'DefaultAxesFontSize', 16)
set(groot, 'DefaultLineLineWidth', 2)
plotWidth = 1200;
plotHeight = 800;

% Load file
[file,path,indx] = uigetfile(".mat",'defname', 'res/');
load([path file]);
%load('res/memEmpty_80steps_amPower_v3/model_test_80step_2s_memEmpty_v3_amp_080W_output.mat');


% some parameters
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




%% Plot

% figure(1)
% plot(t,Xsol(1,:))
% xlabel('Time [s]');
% ylabel([labelStates{1}]);
% title([listStates{1}]);
% set(gcf, 'Position',  [100, 100, plotWidth, plotHeight])
% grid on
markerList = ['o','*', '+','x','s','d','h'];
   
% Mode 
figure(1)
hold on
idData =1;
maxNbMode = 0;
legendTxt = cell(nbControlItem,1);
for i =1:PAR.arch.algos.nb
    legendTxt{i} = ['Algo ' num2str(i) ];
    nbMode = PAR.arch.algos.nbMode(i);
    maxNbMode = max([maxNbMode, nbMode]);
    mode = reshape(result.x(idData:idData+N*nbMode-1),[nbMode,N])'*[1:nbMode]';
    idData = idData + N*nbMode;
    plot(mode,['-' markerList(i)], 'MarkerIndices', 5:5:N)
end

for i =1:PAR.arch.sensors.nb
    legendTxt{PAR.arch.algos.nb+i} = ['Sensor ' num2str(i)];
    nbMode = PAR.arch.sensors.nbMode(i);
    maxNbMode = max([maxNbMode, nbMode]);
    mode = reshape(result.x(idData:idData +N*nbMode-1),[nbMode,N])'*[1:nbMode]';
    idData = idData + N*nbMode;
    plot(mode,['-' markerList(PAR.arch.algos.nb+i)], 'MarkerIndices', 5:5:N)
end

legendTxt{end} = 'OBC';
nbMode = PAR.arch.OBC.nbMode;
maxNbMode = max([maxNbMode, nbMode]);
mode = reshape(result.x(idData:idData +N*nbMode-1),[nbMode,N])'*[1:nbMode]';
plot(mode,['-' markerList(end)], 'MarkerIndices', 5:5:N)


title("Mode")
set(gcf, 'Position',  [100, 100, plotWidth, plotHeight])
ylabel('Actual Mode');
xlabel('Step');
ylim([0.5;maxNbMode+0.5]);
set(gca,'YTick',[0:maxNbMode+1])
grid on
legend(legendTxt)
hold off

%Mode multi window

figure(2)
idData =1;
maxNbMode = 0;
legendTxt = cell(PAR.arch.algos.nb,1);
subplot(3,1,1)
title("Algorithm mode")
hold on
for i =1:PAR.arch.algos.nb
    legendTxt{i} = ['Algo ' num2str(i) ];
    nbMode = PAR.arch.algos.nbMode(i);
    maxNbMode = max([maxNbMode, nbMode]);
    mode = reshape(result.x(idData:idData+N*nbMode-1),[nbMode,N])'*[1:nbMode]';
    idData = idData + N*nbMode;
    plot(mode,['-' markerList(i)], 'MarkerIndices', 5:5:N)
end
hold off
ylabel('Mode');
xlabel('Step');
ylim([0.5;maxNbMode+0.5]);
set(gca,'YTick',[0:maxNbMode+1])
grid on
legend(legendTxt)


subplot(3,1,2)
legendTxt = cell(PAR.arch.sensors.nb,1);
hold on
title("Sensor mode")
for i =1:PAR.arch.sensors.nb
    legendTxt{i} = ['Sensor ' num2str(i)];
    nbMode = PAR.arch.sensors.nbMode(i);
    maxNbMode = max([maxNbMode, nbMode]);
    mode = reshape(result.x(idData:idData +N*nbMode-1),[nbMode,N])'*[1:nbMode]';
    idData = idData + N*nbMode;
    plot(mode,['-' markerList(PAR.arch.algos.nb+i)], 'MarkerIndices', 5:5:N)
end
ylabel('Mode');
xlabel('Step');
ylim([0.5;maxNbMode+0.5]);
set(gca,'YTick',[0:maxNbMode+1])
grid on
legend(legendTxt)
hold off

subplot(3,1,3)
title("OBC Mode")
hold on
nbMode = PAR.arch.OBC.nbMode;
maxNbMode = nbMode;
mode = reshape(result.x(idData:idData +N*nbMode-1),[nbMode,N])'*[1:nbMode]';
plot(mode,['-' markerList(end)], 'MarkerIndices', 5:5:N)
ylabel('Mode');
xlabel('Step');
ylim([0.5;maxNbMode+0.5]);
set(gca,'YTick',[0:maxNbMode+1])
grid on

set(gcf, 'Position',  [100, 100, plotWidth, plotHeight])
hold off


% Accuracy
figure(3)
idData = nbCtrlVariable+1;
legendTxt = cell(PAR.arch.typeDataOut.nb,1);
hold on
for i =1:PAR.arch.typeDataOut.nb
    legendTxt{i} = ['Data Type ' num2str(i) ];
    accData = reshape(result.x(idData:idData+N-1),[1,N]);
    idData = idData + N;
    plot(accData,['-' markerList(i)], 'MarkerIndices', 5:5:N)
end
title("Accuracy")
set(gcf, 'Position',  [100, 100, plotWidth, plotHeight])
ylabel('Accuracy [0;1]');
xlabel('Step');
ylim([0;1]);
grid on
legend(legendTxt)
hold off

% Power
figure(4)
legendTxt = {'Power available', 'Power consume Sensors', 'Power consume OBC','Power left'};
powerAvai = PAR.powerGen.offset.value+PAR.powerGen.amplitude.value*cos(2*pi*k/PAR.powerGen.period.value);

powerConsumeSensor = zeros(N,1);
idData =sum(N*PAR.arch.algos.nbMode) + 1;

for i =1:PAR.arch.sensors.nb

    nbMode = PAR.arch.sensors.nbMode(i);
    powerConsumeSensor = powerConsumeSensor+ reshape(result.x(idData:idData +N*nbMode-1),[nbMode,N])'*PAR.arch.sensors.power.value{i}';

    idData = idData + N*nbMode;
end

nbMode = PAR.arch.OBC.nbMode;
powerConsumeOBC = reshape(result.x(idData:idData +N*nbMode-1),[nbMode,N])'*PAR.arch.OBC.power.value;

hold on
plot(powerAvai,['-' markerList(1)], 'MarkerIndices', 5:5:N)
plot(powerConsumeSensor,['-' markerList(2)], 'MarkerIndices', 5:5:N)
plot(powerConsumeOBC,['-' markerList(3)], 'MarkerIndices', 5:5:N)
plot(powerAvai-(powerConsumeOBC+powerConsumeSensor)',['-' markerList(4)], 'MarkerIndices', 1:10:N)
title("Power")
set(gcf, 'Position',  [100, 100, plotWidth, plotHeight])
ylabel('Power [W]');
xlabel('Step');
ylim([0;(PAR.powerGen.offset.value+PAR.powerGen.amplitude.value)*1.1]);
grid on
legend(legendTxt)
hold off


% Processing Res
figure(5)
legendTxt = {'Processing resources available', 'Processing resources consumed', 'Processing resources left'};
procRessConsume = zeros(N,1);

idData =1;
for i =1:PAR.arch.algos.nb

    nbMode = PAR.arch.algos.nbMode(i);
    procRessConsume = procRessConsume+ reshape(result.x(idData:idData+N*nbMode-1),[nbMode,N])'*PAR.arch.algos.procRess.value{i}';
    idData = idData + N*nbMode;
end


idData =sum(N*[PAR.arch.algos.nbMode  PAR.arch.sensors.nbMode]) + 1;
nbMode = PAR.arch.OBC.nbMode;
procRessAvai = reshape(result.x(idData:idData +N*nbMode-1),[nbMode,N])'*PAR.arch.OBC.procRess.value;

hold on
plot(procRessAvai,['-' markerList(1)], 'MarkerIndices', 5:5:N)
plot(procRessConsume,['-' markerList(2)], 'MarkerIndices', 5:5:N)
plot(procRessAvai-procRessConsume,['-' markerList(3)], 'MarkerIndices', 5:5:N)
title("Processing resources")
set(gcf, 'Position',  [100, 100, plotWidth, plotHeight])
ylabel('Processing resources [MIPS]');
xlabel('Step');
ylim([0;max(PAR.arch.OBC.procRess.value)*1.1]);
grid on
legend(legendTxt)
hold off


% Memory storage
figure(6)
hold on
%subplot(2,1,1)
idData = sum(sizeVar(1:nbControlItem+4))+1;
legendTxt = cell(PAR.arch.memories.nb,1);
hold on
for i =1:PAR.arch.memories.nb
    legendTxt{i} = ['Memory ' num2str(i) ];
    memUsage = reshape(result.x(idData:idData+N-1),[1,N]);
    idData = idData + N;
    plot(memUsage,['-' markerList(i)], 'MarkerIndices', 5:5:N)
end
title("Memory usage")
set(gcf, 'Position',  [100, 100, plotWidth, plotHeight])
ylabel('Memeory usage [Mb]');
xlabel('Step');
%ylim([0;max(PAR.arch.memories.maxStorage.value)]);
grid on
legend(legendTxt);
legendTxtMemStorage = legendTxt;
ax1 = gca;



% Data rate @ mem
%subplot(2,1,2)
figure(7)
hold on
legendTxt = cell(PAR.arch.memories.nb*2,1);

memeDataRate = zeros(N,PAR.arch.memories.nb);
idData =sum(N*PAR.arch.algos.nbMode) + 1;

for i =1:PAR.arch.sensors.nb
    nbMode = PAR.arch.sensors.nbMode(i);
    memLink = PAR.arch.sensors.memoryLink(i);
    
    memeDataRate(:,memLink) = memeDataRate(:,memLink)+ reshape(result.x(idData:idData +N*nbMode-1),[nbMode,N])'*PAR.arch.sensors.dataRate.value{i}';

    idData = idData + N*nbMode;
end


for i = 1:PAR.arch.memories.nb
    plot(memeDataRate(:,i),['-' markerList((i-1)*2+1)], 'MarkerIndices', 5:5:N)
    plot(ones(N,1).*PAR.arch.memories.maxDataRate.value(i),['-' markerList((i-1)*2+2)], 'MarkerIndices', 5:5:N)
    legendTxt{(i-1)*2+1} = ['Mem ' num2str(i)  ' data rate'];
    legendTxt{(i-1)*2+2} = ['Mem ' num2str(i)  ' max data rate'];    
end

title("Data rate at Memory")
set(gcf, 'Position',  [100, 100, plotWidth, plotHeight])
ylabel('Data rate [Mb/s]');
xlabel('Step');
ylim([0;max(PAR.arch.memories.maxDataRate.value)*1.1]);
grid on
lgd2 = legend(legendTxt);
legendTxtMemRate = legendTxt;
hold off
ax2 = gca;

fnew = figure; % https://ch.mathworks.com/help/matlab/ref/subplot.html
ax1_copy = copyobj(ax1,fnew);
subplot(2,1,1,ax1_copy)
legend(legendTxtMemStorage)

copies = copyobj(ax2,fnew);
ax2_copy = copies(1);
subplot(2,1,2,ax2_copy)
legend(legendTxtMemRate)
set(gcf, 'Position',  [100, 100, plotWidth, plotHeight])


%% ----------------- debug --------------------------------

plotDebugGraph =0;
if plotDebugGraph==1

    idData = sum(sizeVar(1:nbControlItem+1))+1;
    figure 
    hold on
    for i =1:PAR.arch.algos.nb
        mode = reshape(result.x(idData:idData+N-1),[1,N]);
        idData = idData + N;
        plot(mode)
    end
    title("Inc")
    legend
    hold off

    figure 
    hold on
    for i =1:PAR.arch.algos.nb
        mode = reshape(result.x(idData:idData+N-1),[1,N]);
        idData = idData + N;
        plot(mode)
    end
    title("LastUp")
    legend
    hold off

    figure 
    hold on
    for i =1:PAR.arch.algos.nb
        mode = reshape(result.x(idData:idData+N-1),[1,N]);
        idData = idData + N;
        plot(mode)
    end
    title("Algo up")
    legend
    hold off

    figure 
    hold on
    for i =1:PAR.arch.memories.nb
        mode = reshape(result.x(idData:idData+N-1),[1,N]);
        idData = idData + N;
        plot(mode)
    end
    title("Mem usage")
    legend
    hold off

    figure 
    hold on
    for i =1:PAR.arch.typeDataOut.nb
        mode = reshape(result.x(idData:idData+N-1),[1,N]);
        idData = idData + N;
        plot(mode)
    end
    title("Acc >= 0")
    legend
    hold off

    figure 
    hold on
    for i =1:PAR.arch.algos.nb
        mode = reshape(result.x(idData:idData+N-1),[1,N]);
        idData = idData + N;
        plot(mode)
    end
    title("ctrl algo up 1")
    legend
    hold off

    figure 
    hold on
    for i =1:PAR.arch.algos.nb
        mode = reshape(result.x(idData:idData+N-1),[1,N]);
        idData = idData + N;
        plot(mode)
    end
    title("ctrl algo up 2")
    legend
    hold off

end
    




%% ------------------ OLD ---------------------------
% 
% listStates = {'Information on target', 'Rate of information', 'Power available', 'Processing ressources available', 'Data rate at POBC', 'Storage capacity usage', 'Last update time', 'Accuracy Max'};
% labelStates ={'Accuracy [0;1]', 'Frequency [Hz]', 'Power left [W]', 'Processing ressources left [MIPS]', 'Rate [Mbits/s]', 'Memory [MBits]', 'Time [s]', 'Accuracy [0;1]'};
% 
% figure(1)
% plot(t,Xsol(1,:))
% xlabel('Time [s]');
% ylabel([labelStates{1}]);
% title([listStates{1}]);
% set(gcf, 'Position',  [100, 100, plotWidth, plotHeight])
% grid on
% 
% % figure(2)
% % plot(t,Xsol(2,:))
% % xlabel('Time [s]');
% % ylabel([labelStates{2}]);
% % title([listStates{2}]);
% % grid on
% 
% figure(3)
% hold on
% for i = 3:8-2
%     subplot(2,2,i-2)
%     plot(t,Xsol(i,:))
%     xlabel('Time [s]');
%     ylabel([labelStates{i}]);
%     title([listStates{i}]);
%     grid on
% end
% set(gcf, 'Position',  [100, 100, plotWidth, plotHeight])
% hold off
% 
% % figure(4)
% % plot(t,Xsol(7,:))
% % xlabel('Time [s]');
% % ylabel([labelStates{7}]);
% % title([listStates{7}]);
% % grid on
% 
% markerList = ['o','*', '+','x','s','d','h'];
% internVar = 0;
% if internVar==1
%     listVarIntern = {'Frequency algo', 'Accuracy Algo', 'Processing ressources algo', 'Data rate algo'};
%     
%     for i= 1:4
%         figure
%         hold on
%         for j=1:PAR.nb_algo
%             plot(t,Xsol(9+(i-1)*PAR.nb_algo+(j-1),:),['-' markerList(j)])
%         end  
%         hold off
%         title([listVarIntern{i}]);
%     end
% end
% 
% figure
% hold on
% subplot(2,2,1)
% p = plot(t(1:end-1), sensorOn);
% for i = 1:PAR.nb_sensor
%     p(i).Marker = markerList(i);
% end
% ylim([0;1]);
% xlabel('Time [s]');
% ylabel('On/Off');
% title("Sensor Status");
% set(gca,'YTick',[0 1])
% legend()
% grid on
% 
% subplot(2,2,2)
% p = plot(t(1:end-1),sensorMode);
% for i = 1:PAR.nb_sensor
%     p(i).Marker = markerList(i);
% end
% ylim([1;PAR.nb_mode_sensor]);
% xlabel('Time [s]');
% ylabel('Mode');
% title("Sensor mode");
% set(gca,'YTick',1:PAR.nb_mode_sensor)
% legend()
% grid on
% 
% subplot(2,2,3)
% p = plot(t(1:end-1),algoMode);
% for i = 1:PAR.nb_algo
%     p(i).Marker = markerList(i);
% end
% ylim([1;PAR.nb_mode_algo]);
% xlabel('Time [s]');
% ylabel('Mode');
% title("Algo mode");
% set(gca,'YTick',1:PAR.nb_mode_algo)
% legend()
% grid on
% 
% subplot(2,2,4)
% plot(t(1:end-1),OBCMode)
% ylim([1;PAR.nb_mode_OBC]);
% xlabel('Time [s]');
% ylabel('Mode');
% title("POBC mode");
% set(gca,'YTick',1:PAR.nb_mode_algo)
% grid on
% set(gcf, 'Position',  [100, 100, plotWidth, plotHeight])
% hold off
% 
