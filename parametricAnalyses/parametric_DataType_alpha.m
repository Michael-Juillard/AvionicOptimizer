%% OPTIMIZER 
% Project       Resources Management for Space Avionics
% Author        Michael Juillard
% Date          Apr 2021
% Version       4


%% initialize
clear ;
close all;
clc;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%Path on lab pc
addpath('/opt/matlab/libMatlab/gurobi912/linux64/matlab')
%addpath('/home/espacemj/Documents/MATLAB/Gurobi/gurobi9.1.1_linux64/gurobi911/linux64/matlab')

%% parameters
PAR = loadParamXML("param/param_mem_empty_v2.xml");

list_alpha = 0.001:0.0005:0.004;

nbIteration = length(list_alpha);
for iteration=1:nbIteration
    disp(['----------------------------------ITERATION ' num2str(iteration) '/'  num2str(nbIteration)  ' ---------------------------------']);
    PAR.arch.typeDataOut.alpha(1) = list_alpha(iteration);
    PAR.arch.typeDataOut.alpha(2) = list_alpha(iteration);
    
    
    
    %% Simulation
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

    % sizeLinCstr     = N*[PAR.arch.algos.nb, PAR.arch.sensors.nb, 1, ... % mode continuity constraint (algos, sensors, obc)
    %                     1, 1, PAR.arch.memories.nb, PAR.arch.memories.nb,... % Power, Processign Res, Data Rate, Storage constraint
    %                     PAR.arch.algos.nb ,PAR.arch.algos.nb ,PAR.arch.algos.nb]; % Update variable lower bound, upper bound, equality
    % nbLineCstr      = sum(sizeLinCstr);
    % 
    % sizeQuadCstr    = [N*PAR.arch.typeDataOut.nb, (N-1)*PAR.arch.typeDataOut.nb, (N-1)*PAR.arch.typeDataOut.nb, ... & accuracy def cstr, accuracy >= 0 cstr 1, accuracy >= 0 cstr 2
    %                    N*PAR.arch.algos.nb, N*PAR.arch.algos.nb , ... % LastUpdate def cstr, Accuracy inc def cstr
    %                    (N-1)*PAR.arch.algos.nb]; % Algo mode chnage cstr
    % nbQuadCstr      = sum(sizeQuadCstr);   

    %% Build model
    model.modelname = ['model_' PAR.nameTest '_' num2str(N) 'step_' num2str(tf) 's_memEmpty_alpha_' num2str(list_alpha(iteration)*10000) 'div10000' ];
    model.modelsense = 'max';

    %% Name variable
    currentIdVar = 0;
    for r= 1:PAR.arch.algos.nb
        for p = 1:N
            for q = 1:PAR.arch.algos.nbMode(r)
                currentIdVar = currentIdVar+1;
                model.varnames{currentIdVar} = sprintf('Algo%d_mode%d_step%d', r,q,p);
            end
        end
    end

    for r= 1:PAR.arch.sensors.nb
        for p = 1:N
            for q = 1:PAR.arch.sensors.nbMode(r)
                currentIdVar = currentIdVar+1;
                model.varnames{currentIdVar} = sprintf('Sensor%d_mode%d_step%d', r,q,p);
            end
        end
    end

    for p = 1:N
        for q = 1:PAR.arch.OBC.nbMode
            currentIdVar = currentIdVar+1;
            model.varnames{currentIdVar} = sprintf('OBC_mode%d_step%d',q,p);
        end
    end


    for q = 1:PAR.arch.typeDataOut.nb
        for p = 1:N
            currentIdVar = currentIdVar+1;
            model.varnames{currentIdVar} = sprintf('Acc_dataType%d_step%d',q,p);
        end
    end

    for q = 1:PAR.arch.algos.nb
        for p = 1:N
            currentIdVar = currentIdVar+1;
            model.varnames{currentIdVar} = sprintf('Inc_acc_Algo%d_step%d',q,p);
        end
    end

    for q = 1:PAR.arch.algos.nb
        for p = 1:N
            currentIdVar = currentIdVar+1;
            model.varnames{currentIdVar} = sprintf('LastUp_Algo%d_step%d',q,p);
        end
    end

    for q = 1:PAR.arch.algos.nb
        for p = 1:N
            currentIdVar = currentIdVar+1;
            model.varnames{currentIdVar} = sprintf('Up_Algo%d_step%d',q,p);
        end
    end

    for q = 1:PAR.arch.memories.nb
        for p = 1:N
            currentIdVar = currentIdVar+1;
            model.varnames{currentIdVar} = sprintf('Storage_memory%d_step%d',q,p);
        end
    end

    for q = 1:PAR.arch.typeDataOut.nb
        for p = 1:N
            currentIdVar = currentIdVar+1;
            model.varnames{currentIdVar} = sprintf('AccGE0_dataType%d_step%d',q,p);
        end
    end

    for q = 1:PAR.arch.algos.nb
        for p = 1:N
            currentIdVar = currentIdVar+1;
            model.varnames{currentIdVar} = sprintf('Var_up_check_lb_Algo%d_step%d',q,p);
        end
    end

    for q = 1:PAR.arch.algos.nb
        for p = 1:N
            currentIdVar = currentIdVar+1;
            model.varnames{currentIdVar} = sprintf('Var_up_check_ub_Algo%d_step%d',q,p);
        end
    end

    %% Variable type
    model.vtype = [repmat('B', nbCtrlVariable, 1); ... % control variable
                   repmat('C', sum(sizeVar(nbControlItem+1:nbControlItem+2)), 1); ... % accuracy & accuracy increment
                   repmat('I', sizeVar(nbControlItem+3), 1); ... %  last update
                   repmat('B', sizeVar(nbControlItem+4), 1); ...% update variable
                   repmat('C', sizeVar(nbControlItem+5), 1); ... %storage variable
                   repmat('B', sum(sizeVar(nbControlItem+6:end)), 1)];  % Accuarcy > 0 var, lower&upper bound for up var

    %% Objective

    %model.obj   = [zeros(nbCtrlVariable, 1);ones(sizeVar(nbControlItem+1),1); zeros(nbVariable-nbCtrlVariable-sizeVar(nbControlItem+1), 1)];
    % Objective: maximize the sum of the square of tha accuracy
    model.obj = zeros(nbVariable, 1);
    model.Q   = sparse(nbVariable, nbVariable);
    % Accuracy alwys max
    for i=1:PAR.arch.typeDataOut.nb
        for w = 1:N
            model.Q(nbCtrlVariable + (i-1)*N + w, nbCtrlVariable + (i-1)*N + w) = PAR.objFuncWeight.value{1}(i);
        end
    end
    %Memory always high
%     for i=1:PAR.arch.memories.nb
%         for w = 1:N
%             model.Q(sum(sizeVar(1:nbControlItem+4)) + (i-1)*N + w, sum(sizeVar(1:nbControlItem+4)) + (i-1)*N + w) = PAR.objFuncWeight.value{2}(i)*1/PAR.arch.memories.maxStorage.value(i);
%         end
%     end

    %% Bound
    model.lb    =  [zeros(sum(sizeVar(1:nbControlItem+2)), 1);...
                    -1*(N+1)*ones(sizeVar(nbControlItem+3), 1);... %For last update variable
                    zeros(sum(sizeVar(nbControlItem+4:end)), 1)];% lower bound
    model.ub    = [ones(nbCtrlVariable, 1); ... % control variable
                   ones(sum(sizeVar(nbControlItem+1:nbControlItem+2)), 1); ... % accuracy & accuracy increment
                   (N+1)*ones(sizeVar(nbControlItem+3), 1); ... % last update
                   ones(sizeVar(nbControlItem+4), 1); ... % update variable
                   reshape((PAR.arch.memories.maxStorage.value'*ones(1,N))', [],1); ... % storage variable
                   ones(sum(sizeVar(nbControlItem+6:end)), 1)]; % Accuarcy >) 0 var, lower&upper bound for up var


    %% Linear Constraint 

    %Matrix A :

    %Continuity constraint
    Mat_cont = zeros(N*nbControlItem,nbVariable);
    nbVarSoFar = 0;
    %Algo
    for id=1:PAR.arch.algos.nb
        i = id;
        nbMode = PAR.arch.algos.nbMode(i);
        range = (id-1)*N + 1:id*N;
        Mat_cont(range,:) = [zeros(N,nbVarSoFar), kron(eye(N),ones(1,nbMode)), zeros(N,nbVariable-nbVarSoFar-N*nbMode)];
        nbVarSoFar = nbVarSoFar+N*nbMode;
    end
    % Sensors
    for id=PAR.arch.algos.nb+1:PAR.arch.algos.nb+PAR.arch.sensors.nb
        i = id-PAR.arch.algos.nb;
        nbMode = PAR.arch.sensors.nbMode(i);
        range = (id-1)*N + 1:id*N;
        Mat_cont(range,:) = [zeros(N,nbVarSoFar), kron(eye(N),ones(1,nbMode)), zeros(N,nbVariable-nbVarSoFar-N*nbMode)];
        nbVarSoFar = nbVarSoFar+N*nbMode;
    end
    %OBC
    nbMode = PAR.arch.OBC.nbMode;
    id = PAR.arch.algos.nb+PAR.arch.sensors.nb+1;
    range = (id-1)*N + 1:id*N;
    Mat_cont(range,:) = [zeros(N,nbVarSoFar), kron(eye(N),ones(1,nbMode)), zeros(N,nbVariable-nbVarSoFar-N*nbMode)];
    Mat_cont= sparse(Mat_cont);

    %Inital mode 
    Mat_mode_init = zeros(nbControlItem,nbVariable);
    elementId = 0;
    for i=1:PAR.arch.algos.nb
        elementId = elementId + 1;
        Mat_mode_init(elementId, sum(sizeVar(1:(elementId-1))) + PAR.arch.algos.modeInit(i) )=1;
    end
    for i=1:PAR.arch.sensors.nb
        elementId = elementId + 1;
        Mat_mode_init(elementId, sum(sizeVar(1:(elementId-1))) + PAR.arch.sensors.modeInit(i) )=1;
    end
    elementId = elementId + 1;
    Mat_mode_init(elementId, sum(sizeVar(1:(elementId-1))) + PAR.arch.OBC.modeInit )=1;
    Mat_mode_init = sparse(Mat_mode_init);


    %Power constraint
    Mat_power = zeros(N,nbVariable);
    nbVarSoFar = sum(sizeVar(1:PAR.arch.algos.nb));
    %sensor
    for i=1:PAR.arch.sensors.nb
        nbMode = PAR.arch.sensors.nbMode(i);
        range = nbVarSoFar + 1:nbVarSoFar+N*nbMode;
        Mat_power(:,range) = kron(eye(N),PAR.arch.sensors.power.value{i});
        nbVarSoFar = nbVarSoFar+N*nbMode;
    end
    %OBC
    nbMode = PAR.arch.OBC.nbMode;
    range = nbVarSoFar + 1:nbVarSoFar+N*nbMode;
    Mat_power(:,range) = kron(eye(N),PAR.arch.OBC.power.value');
    Mat_power= sparse(Mat_power);

    %Processing resources constraint
    Mat_procRes = zeros(N,nbVariable);
    nbVarSoFar = 0;
    %algo
    for i=1:PAR.arch.algos.nb
        nbMode = PAR.arch.algos.nbMode(i);
        range = nbVarSoFar + 1:nbVarSoFar+N*nbMode;
        Mat_procRes(:,range) = -1*kron(eye(N),PAR.arch.algos.procRess.value{i});
        nbVarSoFar = nbVarSoFar+N*nbMode;
    end
    nbVarSoFar = sum(sizeVar(1:(PAR.arch.algos.nb+PAR.arch.sensors.nb)));
    %OBC
    nbMode = PAR.arch.OBC.nbMode;
    range = nbVarSoFar + 1:nbVarSoFar+N*nbMode;
    Mat_procRes(:,range) = kron(eye(N),PAR.arch.OBC.procRess.value');
    Mat_procRes= sparse(Mat_procRes);

    %Data rate constraint (per memory !)
    Mat_dataRate = zeros(N*PAR.arch.memories.nb,nbVariable);
    for j=1:PAR.arch.memories.nb
        nbVarSoFar = sum(sizeVar(1:PAR.arch.algos.nb));
        %sensor
        for i=1:PAR.arch.sensors.nb
            nbMode = PAR.arch.sensors.nbMode(i);
            if PAR.arch.sensors.memoryLink(i)==j %if sensor link to memory
                rangeX = (j-1)*N+1:j*N;
                rangeY = nbVarSoFar + 1:nbVarSoFar+N*nbMode;
                Mat_dataRate(rangeX,rangeY) = kron(eye(N),PAR.arch.sensors.dataRate.value{i});
            end
            nbVarSoFar = nbVarSoFar+N*nbMode;
        end
    end
    Mat_dataRate = sparse(Mat_dataRate);

    %Storage constraint (for every memory)
    Mat_storage = zeros(N*PAR.arch.memories.nb,nbVariable);
    subMat_Storage = eye(N);
    subMat_Storage(2:end,1:end-1) =  subMat_Storage(2:end,1:end-1)-1*eye(N-1);
    for j=1:PAR.arch.memories.nb

        rangeX = (j-1)*N+1:j*N;
        %Algo
        nbVarAlgoSoFar = 0;
        for i=1:PAR.arch.algos.nb
            nbMode = PAR.arch.algos.nbMode(i);
            if PAR.arch.algos.memoryLink(i)==j %if algo link to memory
                rangeY = nbVarAlgoSoFar + 1:nbVarAlgoSoFar+N*nbMode;
                Mat_storage(rangeX,rangeY) = dt*kron(eye(N),PAR.arch.algos.procDataRate.value{i});
            end
            nbVarAlgoSoFar = nbVarAlgoSoFar+N*nbMode;
        end

        %sensor
        nbVarSensorSoFar = nbVarAlgoSoFar;
        for i=1:PAR.arch.sensors.nb
            nbMode = PAR.arch.sensors.nbMode(i);
            if PAR.arch.sensors.memoryLink(i)==j %if sensor link to memory
                rangeY = nbVarSensorSoFar + 1:nbVarSensorSoFar+N*nbMode;
                Mat_storage(rangeX,rangeY) = -dt*kron(eye(N),PAR.arch.sensors.dataRate.value{i});
            end
            nbVarSensorSoFar = nbVarSensorSoFar+N*nbMode;
        end

        %Storage variable
        rangeY = sum(sizeVar(1:nbControlItem+4))+(j-1)*N+1:sum(sizeVar(1:nbControlItem+4))+j*N;
        Mat_storage(rangeX,rangeY) = subMat_Storage;

        %replacement of 1st condition (only Storage variable)
        Mat_storage((j-1)*N+1,1:nbCtrlVariable) = zeros(1,nbCtrlVariable);
    end
    Mat_storage = sparse(Mat_storage);

    %Lower bound of the update variable constraint (for every algo) no step 1
    Mat_lb_UpVar = zeros((N-1)*PAR.arch.algos.nb,nbVariable);
    nbVarAlgoSoFar = 0;
    for i=1:PAR.arch.algos.nb
        rangeX = (i-1)*(N-1)+1:i*(N-1);

        %Algo
        nbMode = PAR.arch.algos.nbMode(i);
        rangeY = nbVarAlgoSoFar + 1:nbVarAlgoSoFar+N*nbMode;
        Mat_lb_UpVar(rangeX,rangeY) = [ zeros(N-1,nbMode), -1*kron(eye(N-1),PAR.arch.algos.period.value{i})];
        nbVarAlgoSoFar = nbVarAlgoSoFar+N*nbMode;

        % lastUpdate (at k-1)
        rangeY = sum(sizeVar(1:nbControlItem+2))+(i-1)*N+1:sum(sizeVar(1:nbControlItem+2))+i*N;
        Mat_lb_UpVar(rangeX,rangeY) = [-1*eye(N-1), zeros(N-1,1)];

    %     % lower boundry var => needed
    %     rangeY = sum(sizeVar(1:nbControlItem+6))+(i-1)*N+1:sum(sizeVar(1:nbControlItem+6))+i*N;
    %     Mat_lb_UpVar(rangeX,rangeY) = [ zeros(N-1,1) PAR.arch.algos.maxPeriod(i)*eye(N-1)];

        % upper boundry var
        rangeY = sum(sizeVar(1:nbControlItem+7))+(i-1)*N+1:sum(sizeVar(1:nbControlItem+7))+i*N;
        Mat_lb_UpVar(rangeX,rangeY) = [ zeros(N-1,1), delta*eye(N-1)];
    end
    Mat_lb_UpVar = sparse(Mat_lb_UpVar);

    %Upper bound of the update variable constraint (for every algo) no step 1
    Mat_ub_UpVar = zeros((N-1)*PAR.arch.algos.nb,nbVariable);
    nbVarAlgoSoFar = 0;
    for i=1:PAR.arch.algos.nb
        rangeX = (i-1)*(N-1)+1:i*(N-1);

        %Algo
        nbMode = PAR.arch.algos.nbMode(i);
        rangeY = nbVarAlgoSoFar + 1:nbVarAlgoSoFar+N*nbMode;
        Mat_ub_UpVar(rangeX,rangeY) =  [ zeros(N-1,nbMode), kron(eye(N-1),PAR.arch.algos.period.value{i})];
        nbVarAlgoSoFar = nbVarAlgoSoFar+N*nbMode;

        % lastUpdate at k-1
        rangeY = sum(sizeVar(1:nbControlItem+2))+(i-1)*N+1:sum(sizeVar(1:nbControlItem+2))+i*N;
        Mat_ub_UpVar(rangeX,rangeY) = [eye(N-1), zeros(N-1,1)];

        % lower boundry var
        rangeY = sum(sizeVar(1:nbControlItem+6))+(i-1)*N+1:sum(sizeVar(1:nbControlItem+6))+i*N;
        Mat_ub_UpVar(rangeX,rangeY) = [ zeros(N-1,1), delta*eye(N-1)];

        % upper boundry var
        rangeY = sum(sizeVar(1:nbControlItem+7))+(i-1)*N+1:sum(sizeVar(1:nbControlItem+7))+i*N;
        Mat_ub_UpVar(rangeX,rangeY) = [ zeros(N-1,1), -1*PAR.arch.algos.maxPeriod(i)*eye(N-1)];
    end
    Mat_ub_UpVar = sparse(Mat_ub_UpVar);

    %Equality of the update variable constraint
    Mat_eq_UpVar = zeros(N*PAR.arch.algos.nb,nbVariable);
    for i=1:PAR.arch.algos.nb
        rangeX = (i-1)*N+1:i*N;

        % up variable
        rangeY = sum(sizeVar(1:nbControlItem+3))+(i-1)*N+1:sum(sizeVar(1:nbControlItem+3))+i*N;
        Mat_eq_UpVar(rangeX,rangeY) = eye(N);

        % lower boundry var
        rangeY = sum(sizeVar(1:nbControlItem+6))+(i-1)*N+1:sum(sizeVar(1:nbControlItem+6))+i*N;
        Mat_eq_UpVar(rangeX,rangeY) = eye(N);

        % upper boundry var
        rangeY = sum(sizeVar(1:nbControlItem+7))+(i-1)*N+1:sum(sizeVar(1:nbControlItem+7))+i*N;
        Mat_eq_UpVar(rangeX,rangeY) = eye(N);

        % for 1st condition up=1
        Mat_eq_UpVar((i-1)*N+1, sum(sizeVar(1:nbControlItem+4)):end ) = zeros(1, sum(sizeVar(nbControlItem+5:end))+1 ); 
    end
    Mat_eq_UpVar = sparse(Mat_eq_UpVar);



    model.A = sparse([Mat_cont;Mat_mode_init;Mat_power;Mat_procRes;Mat_dataRate;Mat_storage;Mat_lb_UpVar;Mat_ub_UpVar;Mat_eq_UpVar]);

    % Vector rhs
    model.rhs = [ones(N*nbControlItem,1);... %continuity 
                 ones(nbControlItem,1);...  % inital mode
                 [PAR.powerGen.offset.value+PAR.powerGen.amplitude.value*cos(2*pi*k/PAR.powerGen.period.value)]';... % power generation
                 zeros(N,1);... %processing ressources
                 reshape((PAR.arch.memories.maxDataRate.value'*ones(1,N))', [],1); ... % data rate
                 reshape([PAR.initCond.storageCap.value; zeros(N-1,PAR.arch.memories.nb)], [],1);... % storage constraint with need initial condition
                 reshape((-1*k_no1)'*ones(1,PAR.arch.algos.nb), [],1);... %lowBound upVar constraint
                 reshape((k_no1)'*ones(1,PAR.arch.algos.nb), [],1);... %upBound upVar constraint
                 ones(N*PAR.arch.algos.nb,1);... % equality upVar constraint
                 ]; 

    % linear constraint sense
    model.sense = [repmat('=', N*nbControlItem, 1);... % continuity 
                   repmat('=', nbControlItem, 1);...% inital mode
                   repmat('<', N, 1);... % power generation
                   repmat('>', N, 1);... % processing ressources
                   repmat('<', N*PAR.arch.memories.nb, 1);... % data rate
                   repmat('=', N*PAR.arch.memories.nb, 1);... % storage constraint
                   repmat('<', (N-1)*PAR.arch.algos.nb, 1);... %lowBound upVar constraint
                   repmat('<', (N-1)*PAR.arch.algos.nb, 1);... %upBound upVar constraint
                   repmat('=', N*PAR.arch.algos.nb, 1);...% equality upVar constraint
                   ]; 


    % linear constraint name
    currentIdCstr = 0;
    for j =1:PAR.arch.algos.nb
        for i =1:N
            currentIdCstr = currentIdCstr+1;
            model.constrnames{currentIdCstr} = sprintf('Continuity_mode_Algo%d_cstr_step%d',j, i);
        end
    end
    for j =1:PAR.arch.sensors.nb
        for i =1:N
            currentIdCstr = currentIdCstr+1;
            model.constrnames{currentIdCstr} = sprintf('Continuity_mode_Sensor%d_cstr_step%d',j, i);
        end
    end
    for i =1:N
        currentIdCstr = currentIdCstr+1;
        model.constrnames{currentIdCstr} = sprintf('Continuity_mode_OBC_cstr_step%d', i);
    end

    for i=1:nbControlItem
        currentIdCstr = currentIdCstr+1;
        model.constrnames{currentIdCstr} = sprintf('Inital_mode_element%d', i);
    end

    for i =1:N
        currentIdCstr = currentIdCstr+1;
        model.constrnames{currentIdCstr} = sprintf('Power_cstr_step%d', i);
    end

    for i =1:N
        currentIdCstr = currentIdCstr+1;
        model.constrnames{currentIdCstr} = sprintf('Processing_Ressources_cstr_step%d', i);
    end

    for j =1:PAR.arch.memories.nb
        for i =1:N
            currentIdCstr = currentIdCstr+1;
            model.constrnames{currentIdCstr} = sprintf('Data_rate_memory%d_cstr_step%d',j, i);
        end
    end

    for j =1:PAR.arch.memories.nb
        for i =1:N
            currentIdCstr = currentIdCstr+1;
            model.constrnames{currentIdCstr} = sprintf('Storage_memory%d_cstr_step%d',j, i);
        end
    end

    for j =1:PAR.arch.algos.nb
        for i =2:N
            currentIdCstr = currentIdCstr+1;
            model.constrnames{currentIdCstr} = sprintf('LowerBound_upVar_Algo%d_cstr_step%d',j, i);
        end
    end

    for j =1:PAR.arch.algos.nb
        for i =2:N
            currentIdCstr = currentIdCstr+1;
            model.constrnames{currentIdCstr} = sprintf('UpperBound_upVar_Algo%d_cstr_step%d',j, i);
        end
    end

    for j =1:PAR.arch.algos.nb
        for i =1:N
            currentIdCstr = currentIdCstr+1;
            model.constrnames{currentIdCstr} = sprintf('Equality_upVar_Algo%d_cstr_step%d',j, i);
        end
    end


    %% Quadratic constraint

    currentIdCstr = 0;
    % Accuracy definition constraint (for every type of data)
    % at step 1=> special case ! (equal to init condition )
    for j= 1:PAR.arch.typeDataOut.nb
        for i = 1:N
            currentIdCstr = currentIdCstr+1;
            if i==1
                model.quadcon(currentIdCstr).Qc = sparse(zeros(nbVariable,nbVariable));
                model.quadcon(currentIdCstr).q  = sparse(zeros(nbVariable,1));
                model.quadcon(currentIdCstr).q(nbCtrlVariable+(j-1)*N+1) = 1;
                model.quadcon(currentIdCstr).rhs = PAR.initCond.infoTarget.value(j); %init cond
            else
            %   define for vector q (linear part)
                accMat = zeros(N,1);
                accMat(i-1) = -1;
                accMat(i) = 1; 

                incMat = zeros(N*PAR.arch.algos.nb,1);
                listIncMatQuad = [];
                for l=1:PAR.arch.algos.nb
                    if PAR.arch.algos.typeDataOut.perAlgo(l) == j % if inc link to this type of data
                        incMat( (l-1)*N+i) = -1;
                        listIncMatQuad = [listIncMatQuad, sum(sizeVar(1:nbControlItem+1))+(l-1)*N+i];  %pos of x_{2,i,k} where i = l here
                    end
                end

                accL0Mat = zeros(N,1);
                accL0Mat(i) = PAR.arch.typeDataOut.alpha(j); % not "-" sign !

                %Quad
                model.quadcon(currentIdCstr).Qrow = listIncMatQuad;
                model.quadcon(currentIdCstr).Qcol = (nbCtrlVariable+(j-1)*N+i-1)*ones(1,length(listIncMatQuad)); %pos of x_{1,b,k-1} multiple time
                model.quadcon(currentIdCstr).Qval = ones(1,length(listIncMatQuad)); %simple multiplication

                model.quadcon(currentIdCstr).q  = sparse([zeros(nbCtrlVariable,1);...
                                                          zeros((j-1)*N,1);accMat;zeros(N*PAR.arch.typeDataOut.nb-j*N,1);... % acc update for the right type of data
                                                          incMat;... % algorithm that increase the accuracy
                                                          zeros(sum(sizeVar(nbControlItem+3:nbControlItem+5)),1);...
                                                          zeros((j-1)*N,1);accL0Mat;zeros(N*PAR.arch.typeDataOut.nb-j*N,1);... % accuracy less 0 for the right type of data
                                                          zeros(sum(sizeVar(nbControlItem+7:end)),1)  ]);
                model.quadcon(currentIdCstr).rhs = 0.0;
            end
            model.quadcon(currentIdCstr).sense = '=';
            model.quadcon(currentIdCstr).name = sprintf('Accuracy%d_quadCstr_step%d',j, i);
        end
    end

    % Accuracy greater than 0 constraint 1 (for every type of data)
    % no k = 1
    for j= 1:PAR.arch.typeDataOut.nb
        for i = 2:N
            currentIdCstr = currentIdCstr+1;

            %   define for vector q (linear part)
            accMat = zeros(N,1);
            accMat(i-1) = 1;

            incMat = zeros(N*PAR.arch.algos.nb,1);
            listIncMatQuad = [];
            for l=1:PAR.arch.algos.nb
                if PAR.arch.algos.typeDataOut.perAlgo(l) == j % if inc link to this type of data
                    incMat( (l-1)*N+i) = 1;
                    listIncMatQuad = [listIncMatQuad, sum(sizeVar(1:nbControlItem+1))+(l-1)*N+i];  %pos of x_{2,i,k} where i = l here
                end
            end

            accL0Mat = zeros(N,1);
            accL0Mat(i) = -1*E; % add the M % remove PAR.arch.typeDataOut.alpha(j)

            %Quad
            model.quadcon(currentIdCstr).Qrow = listIncMatQuad;
            model.quadcon(currentIdCstr).Qcol = (nbCtrlVariable+(j-1)*N+i-1)*ones(1,length(listIncMatQuad)); %pos of x_{1,b,k-1} multiple time
            model.quadcon(currentIdCstr).Qval = -1*ones(1,length(listIncMatQuad)); %simple multiplication

            model.quadcon(currentIdCstr).q  = sparse([zeros(nbCtrlVariable,1);...
                                                      zeros((j-1)*N,1);accMat;zeros(N*PAR.arch.typeDataOut.nb-j*N,1);... % acc update for the right type of data
                                                      incMat;... % algorithm that increase the accuracy
                                                      zeros(sum(sizeVar(nbControlItem+3:nbControlItem+5)),1);...
                                                      zeros((j-1)*N,1);accL0Mat;zeros(N*PAR.arch.typeDataOut.nb-j*N,1);... % accuracy less 0 for the right type of data
                                                      zeros(sum(sizeVar(nbControlItem+7:end)),1)  ]);
            model.quadcon(currentIdCstr).rhs = -1*E+PAR.arch.typeDataOut.alpha(j)+delta; 

            model.quadcon(currentIdCstr).sense = '>';
            model.quadcon(currentIdCstr).name = sprintf('Accuracy%d_GE0_quadCstr1_step%d',j, i);
        end
    end

    % Accuracy greater than 0 constraint 2 (for every type of data)
    % no k = 1
    for j= 1:PAR.arch.typeDataOut.nb
        for i = 2:N
            currentIdCstr = currentIdCstr+1;

            %   define for vector q (linear part)
            accMat = zeros(N,1);
            accMat(i-1) = 1;

            incMat = zeros(N*PAR.arch.algos.nb,1);
            listIncMatQuad = [];
            for l=1:PAR.arch.algos.nb
                if PAR.arch.algos.typeDataOut.perAlgo(l) == j % if inc link to this type of data
                    incMat( (l-1)*N+i) = 1;
                    listIncMatQuad = [listIncMatQuad, sum(sizeVar(1:nbControlItem+1))+(l-1)*N+i];  %pos of x_{2,i,k} where i = l here
                end
            end

            accL0Mat = zeros(N,1);
            accL0Mat(i) = -1*E; % add the M % remove PAR.arch.typeDataOut.alpha(j)

            %Quad
            model.quadcon(currentIdCstr).Qrow = listIncMatQuad;
            model.quadcon(currentIdCstr).Qcol = (nbCtrlVariable+(j-1)*N+i-1)*ones(1,length(listIncMatQuad)); %pos of x_{1,b,k-1} multiple time
            model.quadcon(currentIdCstr).Qval = -1*ones(1,length(listIncMatQuad)); %simple multiplication

            model.quadcon(currentIdCstr).q  = sparse([zeros(nbCtrlVariable,1);...
                                                      zeros((j-1)*N,1);accMat;zeros(N*PAR.arch.typeDataOut.nb-j*N,1);... % acc update for the right type of data
                                                      incMat;... % algorithm that increase the accuracy
                                                      zeros(sum(sizeVar(nbControlItem+3:nbControlItem+5)),1);...
                                                      zeros((j-1)*N,1);accL0Mat;zeros(N*PAR.arch.typeDataOut.nb-j*N,1);... % accuracy less 0 for the right type of data
                                                      zeros(sum(sizeVar(nbControlItem+7:end)),1)  ]);
            model.quadcon(currentIdCstr).rhs = PAR.arch.typeDataOut.alpha(j); % smaller than 0

            model.quadcon(currentIdCstr).sense = '<';
            model.quadcon(currentIdCstr).name = sprintf('Accuracy%d_GE0_quadCstr2_step%d',j, i);
        end
    end

    %Last update definition constraint
    % k=1 => lastUp =-1
    for j= 1:PAR.arch.algos.nb
        for i = 1:N
            currentIdCstr = currentIdCstr+1;
            if i==1
                model.quadcon(currentIdCstr).Qc = sparse(zeros(nbVariable,nbVariable));
                model.quadcon(currentIdCstr).q  = sparse(zeros(nbVariable,1));
                model.quadcon(currentIdCstr).q(sum(sizeVar(1:nbControlItem+2))+(j-1)*N+1) = 1;
                model.quadcon(currentIdCstr).rhs = 1;%-1.0*PAR.arch.algos.period.value{j}(PAR.arch.algos.modeInit(j)); % 
            else
                %   define for vector q (linear part)
                lastUpMat = zeros(N,1);
                lastUpMat(i-1) = -1;
                lastUpMat(i) = 1; 

                upMat =  zeros(N,1);
                upMat(i) = -i;

                %Quad
                model.quadcon(currentIdCstr).Qrow = [sum(sizeVar(1:nbControlItem+2))+(j-1)*N+i-1]; %pos of x_{3,i,k-1} 
                model.quadcon(currentIdCstr).Qcol = [sum(sizeVar(1:nbControlItem+3))+(j-1)*N+i]; %pos of x_{4,i,k} 
                model.quadcon(currentIdCstr).Qval = [1]; %simple multiplication

                model.quadcon(currentIdCstr).q  = sparse([zeros(nbCtrlVariable,1);...
                                                          zeros(sizeVar(nbControlItem+1),1);... % acc  
                                                          zeros(sizeVar(nbControlItem+2),1);... % inc
                                                          zeros((j-1)*N,1);lastUpMat;zeros(N*PAR.arch.algos.nb-j*N,1);... %lastUp for the right algo
                                                          zeros((j-1)*N,1);upMat;zeros(N*PAR.arch.algos.nb-j*N,1);... %up var for the right algo
                                                          zeros(sum(sizeVar(nbControlItem+5:end)),1)  ]);
                model.quadcon(currentIdCstr).rhs = 0.0;

            end


            model.quadcon(currentIdCstr).sense = '=';
            model.quadcon(currentIdCstr).name = sprintf('LastUp_def_algo%d_quadCstr_step%d',j, i);
        end
    end

    %Accuracy increment definition constraint
    for j= 1:PAR.arch.algos.nb
        for i = 1:N
            currentIdCstr = currentIdCstr+1;

            %   define for vector q (linear part)
            nbMode = PAR.arch.algos.nbMode(j);
            lastUpMat = zeros(N,1);
            lastUpMat(i) = 1; 

            %Quad
            model.quadcon(currentIdCstr).Qrow = sum(sizeVar(1:j))-nbMode*N +nbMode*(i-1)+[1:nbMode];  %pos of U_{1,i,j,k} for all j 
            model.quadcon(currentIdCstr).Qcol = (sum(sizeVar(1:nbControlItem+3))+(j-1)*N+i)*ones(1,nbMode); %pos of x_{4,i,k} per mode
            model.quadcon(currentIdCstr).Qval = -1*PAR.arch.algos.accuracy.value{j}; %multiply by accuracy

            model.quadcon(currentIdCstr).q  = sparse([zeros(nbCtrlVariable,1);...
                                                      zeros(sizeVar(nbControlItem+1),1);... % acc  
                                                      zeros((j-1)*N,1);lastUpMat;zeros(N*PAR.arch.algos.nb-j*N,1);... % inc
                                                      zeros(sum(sizeVar(nbControlItem+3:end)),1)  ]);
            model.quadcon(currentIdCstr).rhs = 0.0;
            model.quadcon(currentIdCstr).sense = '=';

            model.quadcon(currentIdCstr).name = sprintf('Acc_inc_def_algo%d_quadCstr_step%d',j, i);
        end
    end

    %Algorithm mode change constraint
    % no step 1 
    for j= 1:PAR.arch.algos.nb
        for i = 2:N
            currentIdCstr = currentIdCstr+1;

            %   define for vector q (linear part)
            nbMode = PAR.arch.algos.nbMode(j);
            modeMat = zeros(N*nbMode,1);
            index1 = [1:nbMode]+ nbMode*(i-2);
            index2 = [1:nbMode]+ nbMode*(i-1);
            modeMat(index1) = -1*[0:(nbMode-1)]; %U_{1,i,j,k-1} for all j mutliply by 0 to -nMode+1 (dec)
            modeMat(index2) = [0:(nbMode-1)]; %U_{1,i,j,k} for all j mutliply by 0 to nMode-1 (inc)

            %Quad
            model.quadcon(currentIdCstr).Qrow = sum(sizeVar(1:j))-nbMode*N +nbMode*(i-2)+[1:nbMode*2];  %pos of U_{1,i,j,k-1} for all j then  U_{1,i,j,k} for all j 
            model.quadcon(currentIdCstr).Qcol = (sum(sizeVar(1:nbControlItem+3))+(j-1)*N+i)*ones(1,nbMode*2); %pos of x_{4,i,k} per mode*2
            model.quadcon(currentIdCstr).Qval = [0:(nbMode-1),0:-1:(1-nbMode)]; %mutliply by corresponding factor 0:nMode-1 then 0:-(nbMode-1)

            model.quadcon(currentIdCstr).q  = sparse([zeros(sum(sizeVar(1:j))-nbMode*N,1); modeMat; zeros(sum(sizeVar(1:PAR.arch.algos.nb))-sum(sizeVar(1:j)),1);...
                                                      zeros(nbCtrlVariable-sum(sizeVar(1:PAR.arch.algos.nb)),1);... % the other ctrl variable
                                                      zeros(sum(sizeVar(nbControlItem+1:end)),1)  ]); % the rest
            model.quadcon(currentIdCstr).rhs = 0.0;
            model.quadcon(currentIdCstr).sense = '=';

            model.quadcon(currentIdCstr).name = sprintf('Algo%d_mode_change_quadCstr_step%d',j, i);
        end
    end

    %% Witre model
    gurobi_write(model, ['res/' model.modelname '.lp']);


    %% Solve model

    params.NonConvex = 2;
    params.outputflag = 1;

    %% Test MIP focus
    % https://www.gurobi.com/documentation/9.1/refman/mip_models.html
    % https://www.gurobi.com/documentation/9.1/examples/params_m.html

    % Set a 2 second time limit
    params.TimeLimit = 5;

    % Now solve the model with different values of MIPFocus

    params.MIPFocus = 0;
    result          = gurobi(model, params);
    bestgap         = result.mipgap;
    bestparams      = params;
    for i = 1:3
        params.MIPFocus = i;
        result          = gurobi(model, params);
        if result.mipgap < bestgap
            bestparams = params;
            bestgap    = result.mipgap;
        end
    end

    % Finally, reset the time limit and Re-solve model to optimality
    
    %bestparams.TimeLimit = Inf;
    bestparams.TimeLimit = 60*60*4;
    %params.MIPFocus = 3; % 

    result = gurobi(model, bestparams);

    %% Analyse result
    disp(result);
    plotData=0;
    if strcmp(result.status, 'OPTIMAL')

        fprintf('The optimal objective is %g\n', result.objval);
    %     fprintf('variable:\n');
    %     for v=1:length(model.varnames)
    %         fprintf('%s %d\n', model.varnames{v}, result.x(v));
    %     end
        %save result  !!!!
        nameFile = ['res/' model.modelname '_output.mat']; 
        save(nameFile, 'result', 'PAR');
    elseif strcmp(result.status, 'TIME_LIMIT')
        fprintf('TIME LIMIT SOL : The optimal objective is %g\n', result.objval);

        nameFile = ['res/' model.modelname '_TM_output.mat']; 
        save(nameFile, 'result', 'PAR');
    else

        if strcmp(result.status, 'INFEASIBLE')
            fprintf('Problem is infeasible.... computing IIS\n');


%             prompt = 'Do you want to compute the IIS? Y/N [Y]: ';
%             str = input(prompt,'s');
%             if isempty(str)
%                 str = 'Y';
%             end
            str = 'Y';

            if str=='Y' % need to save IIS !
                iis = gurobi_iis(model, params);

                nameFile = ['res/IIS_' model.modelname '.mat']; 
                save(nameFile, 'iis', 'PAR');

                if iis.minimal
                    fprintf('IIS is minimal\n');
                else
                    fprintf('IIS is not minimal\n');
                end

                if any(iis.Arows)
                    fprintf('Rows in IIS: \n');
                    disp(strjoin(model.constrnames(iis.Arows),' \n'));
                end
                if any(iis.lb)
                    fprintf('LB in IIS: \n');
                    disp(strjoin(model.varnames(iis.lb),' \n'));
                end
                if any(iis.ub)
                    fprintf('UB in IIS: \n');
                    disp(strjoin(model.varnames(iis.ub),' \n'));
                end
            end
        else
            % Just to handle user interruptions or other problems
            fprintf('Unexpected status %s\n',result.status);
        end
    end

end