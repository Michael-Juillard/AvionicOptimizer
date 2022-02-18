function [PARAM] = loadParamXML(fileInfo)
%loadParamXML Read xml file to laod sim parameters
%   test case : PAR = loadParamXML("param/param_test.xml")
%   Reading and extracting parameters from the xml document
    DOMnode = xmlread(fileInfo);
    xRoot = DOMnode.getDocumentElement;
    if xRoot.getTagName ~= 'param'
        return;
    end

    % nameTest
    nameChild = xRoot.getElementsByTagName('name');
    if nameChild.getLength ~= 0
        PARAM.nameTest = char(nameChild.item(0).getTextContent);
    end

    % satellite
    satelliteChild = xRoot.getElementsByTagName('satellite');
    if satelliteChild.getLength ~= 0
        PARAM.satellite = char(satelliteChild.item(0).getTextContent);
    end

    %Architecture 
    architectureChild = xRoot.getElementsByTagName('architecture');
    if architectureChild.getLength ~= 0
        architectureChild = architectureChild.item(0);
        
        %Name arch
        archNameChild = architectureChild.getElementsByTagName('name');
        if archNameChild.getLength ~= 0
            PARAM.arch.name = char(archNameChild.item(0).getTextContent);
        end
        
        %Memories
        memeoriesChild = architectureChild.getElementsByTagName('memories');
        if memeoriesChild.getLength ~= 0
            oneChild = memeoriesChild.item(0);              
            nbMemoriesChild = oneChild.getLength;
            PARAM.arch.memories.nb = str2double(oneChild.getAttribute('nb'));
            
            for i = 1:nbMemoriesChild-1
                oneMem=oneChild.item(i);
                if oneMem.hasAttributes()
                    idMem = str2double(oneMem.getAttribute('id'));

                    %maxDataRate
                    maxDataRateChild = oneMem.getElementsByTagName('maxDataRate');
                    if maxDataRateChild.getLength ~= 0
                        PARAM.arch.memories.maxDataRate.unit(idMem) = {char(maxDataRateChild.item(0).getAttribute('unit'))};
                        PARAM.arch.memories.maxDataRate.value(idMem) = str2double(maxDataRateChild.item(0).getTextContent);
                    end

                    %maxStorage
                    maxStorageChild = oneMem.getElementsByTagName('maxStorage');
                    if maxStorageChild.getLength ~= 0
                        PARAM.arch.memories.maxStorage.unit(idMem) = {char(maxStorageChild.item(0).getAttribute('unit'))};
                        PARAM.arch.memories.maxStorage.value(idMem) = str2double(maxStorageChild.item(0).getTextContent);
                    end   
                end
            end
        end

        %Sensors
        sensorsChild = architectureChild.getElementsByTagName('sensors');
        if sensorsChild.getLength ~= 0
            oneChild = sensorsChild.item(0);              
            nbSensorChild = oneChild.getLength;
            PARAM.arch.sensors.nb = str2double(oneChild.getAttribute('nb'));

            for i = 1:nbSensorChild-1
                oneSensor=oneChild.item(i);
                if oneSensor.hasAttributes()
                    idSensor = str2double(oneSensor.getAttribute('id'));
                    PARAM.arch.sensors.memoryLink(idSensor) = str2double(oneSensor.getAttribute('memLink'));

                    %name sensor
                    nameSensorChild = oneSensor.getElementsByTagName('name');
                    if nameSensorChild.getLength ~= 0
                        PARAM.arch.sensors.name(idSensor) = {char(nameSensorChild.item(0).getTextContent)};
                    end
                    %mode sensor
                    modesChild = oneSensor.getElementsByTagName('modes');
                    if modesChild.getLength ~= 0
                        oneModeChild = modesChild.item(0); 
                        nbModeChild = oneModeChild.getLength;
                        nbMode = str2double(oneModeChild.getAttribute('nbMode'));
                        PARAM.arch.sensors.modeInit(idSensor) = str2double(oneModeChild.getAttribute('modeStart'));
                        PARAM.arch.sensors.nbMode(idSensor) = nbMode;
                        for j = 1:nbModeChild-1
                            oneMode=oneModeChild.item(j);
                            if oneMode.hasAttributes()
                                idMode = str2double(oneMode.getAttribute('id'));
                                %name Mode
                                nameModeChild = oneMode.getElementsByTagName('name');
                                if nameModeChild.getLength ~= 0
                                    PARAM.arch.sensors.nameMode{idSensor}(idMode) = {char(nameModeChild.item(0).getTextContent)};
                                end
                                % power
                                powerModeChild = oneMode.getElementsByTagName('power');
                                if powerModeChild.getLength ~= 0
                                     PARAM.arch.sensors.power.unit{idSensor}(idMode) = {char(powerModeChild.item(0).getAttribute('unit'))};
                                     PARAM.arch.sensors.power.value{idSensor}(idMode) = str2double(powerModeChild.item(0).getTextContent);
                                end
                                % dataRate
                                dataRateChild = oneMode.getElementsByTagName('dataRate');
                                if dataRateChild.getLength ~= 0
                                     PARAM.arch.sensors.dataRate.unit{idSensor}(idMode) = {char(dataRateChild.item(0).getAttribute('unit'))};
                                     PARAM.arch.sensors.dataRate.value{idSensor}(idMode) = str2double(dataRateChild.item(0).getTextContent);
                                end
                            end
                        end
                    end
                end
            end
        end
        
        %Type Data out
        typeDataOutChild = architectureChild.getElementsByTagName('typeDataOut');
        if typeDataOutChild.getLength ~= 0
            oneChild = typeDataOutChild.item(0);              
            nbTypeDataOutChild = oneChild.getLength;
            PARAM.arch.typeDataOut.nb = str2double(oneChild.getAttribute('nb'));
            
            for i = 1:nbTypeDataOutChild-1
                oneType=oneChild.item(i);
                if oneType.hasAttributes()
                    idType = str2double(oneType.getAttribute('id'));
                    
                    %name dataOut
                    nameDataOutChild = oneType.getElementsByTagName('name');
                    if nameDataOutChild.getLength ~= 0
                        PARAM.arch.typeDataOut.name(idType) = {char(nameDataOutChild.item(0).getTextContent)};
                    end
                    
                    %alpha
                    alphaChild = oneType.getElementsByTagName('alpha');
                    if alphaChild.getLength ~= 0
                        PARAM.arch.typeDataOut.alpha(idType) = str2double(alphaChild.item(0).getTextContent);
                    end

                end
            end
        end
        
        %Algo
        algosChild = architectureChild.getElementsByTagName('algos');
        if algosChild.getLength ~= 0
            oneChild = algosChild.item(0);              
            nbAlgoChild = oneChild.getLength;
            PARAM.arch.algos.nb = str2double(oneChild.getAttribute('nb'));

            for i = 1:nbAlgoChild-1
                oneAlgo=oneChild.item(i);
                if oneAlgo.hasAttributes()
                    idAlgo = str2double(oneAlgo.getAttribute('id'));
                    PARAM.arch.algos.memoryLink(idAlgo) = str2double(oneAlgo.getAttribute('memLink'));
                    PARAM.arch.algos.typeDataOut.perAlgo(idAlgo) = str2double(oneAlgo.getAttribute('typeData'));

                    %name algo
                    nameAlgoChild = oneAlgo.getElementsByTagName('name');
                    if nameAlgoChild.getLength ~= 0
                        PARAM.arch.algos.name(idAlgo) = {char(nameAlgoChild.item(0).getTextContent)};
                    end
                    %mode algo
                    modesChild = oneAlgo.getElementsByTagName('modes');
                    if modesChild.getLength ~= 0
                        oneModeChild = modesChild.item(0); 
                        nbModeChild = oneModeChild.getLength;
                        nbMode = str2double(oneModeChild.getAttribute('nbMode'));
                        PARAM.arch.algos.modeInit(idAlgo) = str2double(oneModeChild.getAttribute('modeStart'));
                        PARAM.arch.algos.nbMode(idAlgo) = nbMode;
                        for j = 1:nbModeChild-1
                            oneMode=oneModeChild.item(j);
                            if oneMode.hasAttributes()
                                idMode = str2double(oneMode.getAttribute('id'));
                                %name Mode
                                nameModeChild = oneMode.getElementsByTagName('name');
                                if nameModeChild.getLength ~= 0
                                    PARAM.arch.algos.nameMode{idAlgo}(idMode) = {char(nameModeChild.item(0).getTextContent)};
                                end
                                % period
                                periodModeChild = oneMode.getElementsByTagName('period');
                                if periodModeChild.getLength ~= 0
                                     PARAM.arch.algos.period.unit{idAlgo}(idMode) = {char(periodModeChild.item(0).getAttribute('unit'))};
                                     PARAM.arch.algos.period.value{idAlgo}(idMode) = str2double(periodModeChild.item(0).getTextContent);
                                end
                                % accuracy
                                accuracyChild = oneMode.getElementsByTagName('accuracy');
                                if accuracyChild.getLength ~= 0
                                     PARAM.arch.algos.accuracy.unit{idAlgo}(idMode) = {char(accuracyChild.item(0).getAttribute('unit'))};
                                     PARAM.arch.algos.accuracy.value{idAlgo}(idMode) = str2double(accuracyChild.item(0).getTextContent);
                                end
                                % procRess
                                procRessChild = oneMode.getElementsByTagName('procRess');
                                if procRessChild.getLength ~= 0
                                     PARAM.arch.algos.procRess.unit{idAlgo}(idMode) = {char(procRessChild.item(0).getAttribute('unit'))};
                                     PARAM.arch.algos.procRess.value{idAlgo}(idMode) = str2double(procRessChild.item(0).getTextContent);
                                end
                                % procDataRate
                                procDataRateChild = oneMode.getElementsByTagName('procDataRate');
                                if procDataRateChild.getLength ~= 0
                                     PARAM.arch.algos.procDataRate.unit{idAlgo}(idMode) = {char(procDataRateChild.item(0).getAttribute('unit'))};
                                     PARAM.arch.algos.procDataRate.value{idAlgo}(idMode) = str2double(procDataRateChild.item(0).getTextContent);
                                end
                            end
                        end
                        PARAM.arch.algos.maxPeriod(idAlgo) = max(PARAM.arch.algos.period.value{idAlgo});
                    end
                end
            end
        end
        
        %OBC
        obcChild = architectureChild.getElementsByTagName('obc');
        if obcChild.getLength ~= 0
            oneOBC = obcChild.item(0);              
 
            %mode OBC
            modesChild = oneOBC.getElementsByTagName('modes');
            if modesChild.getLength ~= 0
                oneModeChild = modesChild.item(0); 
                nbModeChild = oneModeChild.getLength;
                nbMode = str2double(oneModeChild.getAttribute('nbMode'));
                PARAM.arch.OBC.modeInit = str2double(oneModeChild.getAttribute('modeStart'));
                PARAM.arch.OBC.nbMode = nbMode;
                for i = 1:nbModeChild-1
                    oneMode=oneModeChild.item(i);
                    if oneMode.hasAttributes()
                        idMode = str2double(oneMode.getAttribute('id'));
                        %name Mode
                        nameModeChild = oneMode.getElementsByTagName('name');
                        if nameModeChild.getLength ~= 0
                            PARAM.arch.OBC.nameMode(idMode,1) = {char(nameModeChild.item(0).getTextContent)};
                        end
                        % power
                        powerModeChild = oneMode.getElementsByTagName('power');
                        if powerModeChild.getLength ~= 0
                             PARAM.arch.OBC.power.unit(idMode,1) = {char(powerModeChild.item(0).getAttribute('unit'))};
                             PARAM.arch.OBC.power.value(idMode,1) = str2double(powerModeChild.item(0).getTextContent);
                        end
                        % procRess
                        procRessChild = oneMode.getElementsByTagName('procRess');
                        if procRessChild.getLength ~= 0
                             PARAM.arch.OBC.procRess.unit(idMode,1) = {char(procRessChild.item(0).getAttribute('unit'))};
                             PARAM.arch.OBC.procRess.value(idMode,1) = str2double(procRessChild.item(0).getTextContent);
                        end
                    end
                end
            end

        end
        
    end  % End of arch
    
    % alphaParam [move to type data out]
%     alphaParamChild = xRoot.getElementsByTagName('alpha');
%     if alphaParamChild.getLength ~= 0
%         PARAM.alpha = str2double(alphaParamChild.item(0).getTextContent);
%     end
    
    %powerGeneration
    powerGenerationChild = xRoot.getElementsByTagName('powerGeneration');
    if powerGenerationChild.getLength ~= 0
        powerGenerationChild = powerGenerationChild.item(0);
        
        %type
        typeChild = powerGenerationChild.getElementsByTagName('type');
        if typeChild.getLength ~= 0
            PARAM.powerGen.type = char(typeChild.item(0).getTextContent);
        end
        
        %offset
        offsetChild = powerGenerationChild.getElementsByTagName('offset');
        if offsetChild.getLength ~= 0
            PARAM.powerGen.offset.unit = char(offsetChild.item(0).getAttribute('unit'));
            PARAM.powerGen.offset.value = str2double(offsetChild.item(0).getTextContent);
        end
        
        %amplitude
        amplitudeChild = powerGenerationChild.getElementsByTagName('amplitude');
        if amplitudeChild.getLength ~= 0
            PARAM.powerGen.amplitude.unit = char(amplitudeChild.item(0).getAttribute('unit'));
            PARAM.powerGen.amplitude.value = str2double(amplitudeChild.item(0).getTextContent);
        end
        
        %period
        periodChild = powerGenerationChild.getElementsByTagName('period');
        if periodChild.getLength ~= 0
            PARAM.powerGen.period.unit = char(periodChild.item(0).getAttribute('unit'));
            PARAM.powerGen.period.value = str2double(periodChild.item(0).getTextContent);
        end
        
    end
    
%     %constraints [move to memories]
%     constraintsChild = xRoot.getElementsByTagName('constraints');
%     if constraintsChild.getLength ~= 0
%         constraintsChild = constraintsChild.item(0);
%                
%         %maxDataRate
%         maxDataRateChild = constraintsChild.getElementsByTagName('maxDataRate');
%         if maxDataRateChild.getLength ~= 0
%             PARAM.constraints.maxDataRate.unit = char(maxDataRateChild.item(0).getAttribute('unit'));
%             PARAM.constraints.maxDataRate.value = str2double(maxDataRateChild.item(0).getTextContent);
%         end
%         
%         %maxStorage
%         maxStorageChild = constraintsChild.getElementsByTagName('maxStorage');
%         if maxStorageChild.getLength ~= 0
%             PARAM.constraints.maxStorage.unit = char(maxStorageChild.item(0).getAttribute('unit'));
%             PARAM.constraints.maxStorage.value = str2double(maxStorageChild.item(0).getTextContent);
%         end       
%     end
    
    %simulationParameters
    simulationParametersChild = xRoot.getElementsByTagName('simulationParameters');
    if simulationParametersChild.getLength ~= 0
        simulationParametersChild = simulationParametersChild.item(0);
               
        %startTime
        startTimeChild = simulationParametersChild.getElementsByTagName('startTime');
        if startTimeChild.getLength ~= 0
            PARAM.sim.startTime.unit = char(startTimeChild.item(0).getAttribute('unit'));
            PARAM.sim.startTime.value = str2double(startTimeChild.item(0).getTextContent);
        end
        
        %endTime
        endTimeChild = simulationParametersChild.getElementsByTagName('endTime');
        if endTimeChild.getLength ~= 0
            PARAM.sim.endTime.unit = char(endTimeChild.item(0).getAttribute('unit'));
            PARAM.sim.endTime.value = str2double(endTimeChild.item(0).getTextContent);
        end       
        
        %step
        stepChild = simulationParametersChild.getElementsByTagName('step');
        if stepChild.getLength ~= 0
            PARAM.sim.step.unit = char(stepChild.item(0).getAttribute('unit'));
            PARAM.sim.step.value = str2double(stepChild.item(0).getTextContent);
        end  
    end
    
    %initalConditionState
    initCondChild = xRoot.getElementsByTagName('initalConditionState');
    if initCondChild.getLength ~= 0
        initCondChild = initCondChild.item(0);
               
        %infoTarget
        infoTargetChild = initCondChild.getElementsByTagName('infoTarget');
        if infoTargetChild.getLength ~= 0
            PARAM.initCond.infoTarget.name = char(infoTargetChild.item(0).getAttribute('name'));
                        
            oneChild = infoTargetChild.item(0);              
            nbInfoTargetChild = oneChild.getLength;           
            for i = 1:nbInfoTargetChild-1
                oneInfoTarget=oneChild.item(i);
                if oneInfoTarget.hasAttributes()
                    id = str2double(oneInfoTarget.getAttribute('id'));
                    PARAM.initCond.infoTarget.value(id) = str2double(oneInfoTarget.item(0).getTextContent);
                end
            end
        end
        
        %power
        freqChild = initCondChild.getElementsByTagName('avgFrequency');
        if freqChild.getLength ~= 0
            PARAM.initCond.avgFrequency.name = char(freqChild.item(0).getAttribute('name'));
            PARAM.initCond.avgFrequency.value = str2double(freqChild.item(0).getTextContent);
        end  
        
        %power
        powerChild = initCondChild.getElementsByTagName('power');
        if powerChild.getLength ~= 0
            PARAM.initCond.power.name = char(powerChild.item(0).getAttribute('name'));
            PARAM.initCond.power.value = str2double(powerChild.item(0).getTextContent);
        end       
        
        %procRess
        procRessChild = initCondChild.getElementsByTagName('procRess');
        if procRessChild.getLength ~= 0
            PARAM.initCond.procRess.name = char(procRessChild.item(0).getAttribute('name'));
            PARAM.initCond.procRess.value = str2double(procRessChild.item(0).getTextContent);
        end  
        
        %dataRate
        dataRateChild = initCondChild.getElementsByTagName('dataRate');
        if dataRateChild.getLength ~= 0
            PARAM.initCond.dataRate.name = char(dataRateChild.item(0).getAttribute('name'));
            
            oneChild = dataRateChild.item(0);              
            nbDataRateChild = oneChild.getLength;           
            for i = 1:nbDataRateChild-1
                oneDataRate=oneChild.item(i);
                if oneDataRate.hasAttributes()
                    id = str2double(oneDataRate.getAttribute('id'));
                    PARAM.initCond.dataRate.value(id) = str2double(oneDataRate.item(0).getTextContent);
                end
            end
        end  
        
        %storageCap
        storageCapChild = initCondChild.getElementsByTagName('storageCap');
        if storageCapChild.getLength ~= 0
            PARAM.initCond.storageCap.name = char(storageCapChild.item(0).getAttribute('name'));
            PARAM.initCond.storageCap.value = str2double(storageCapChild.item(0).getTextContent);
            
            oneChild = storageCapChild.item(0);              
            nbStorageCapChild = oneChild.getLength;           
            for i = 1:nbStorageCapChild-1
                oneStorageCap=oneChild.item(i);
                if oneStorageCap.hasAttributes()
                    id = str2double(oneStorageCap.getAttribute('id'));
                    PARAM.initCond.storageCap.value(id) = str2double(oneStorageCap.item(0).getTextContent);
                end
            end
        end  
    end
    
    
    %objectiveFunctionWeights
    objectiveFunctionWeightsChild = xRoot.getElementsByTagName('objectiveFunctionWeights');
    if objectiveFunctionWeightsChild.getLength ~= 0
        oneWeightChild = objectiveFunctionWeightsChild.item(0);   
        nbWeightChild = oneWeightChild.getLength;
        for i = 1:nbWeightChild-1
            oneWeight=oneWeightChild.item(i);
            if oneWeight.hasAttributes()
                idWeight = str2double(oneWeight.getAttribute('id'));
                PARAM.objFuncWeight.name(idWeight) = {char(oneWeight.getAttribute('name'))};
                
                nbValueChild = oneWeight.getLength;
                for j = 1:nbValueChild-1
                    oneValue = oneWeight.item(j);
                    if oneValue.hasAttributes()
                        idValue = str2double(oneValue.getAttribute('id'));
                        PARAM.objFuncWeight.value{idWeight}(idValue) = str2double(oneValue.item(0).getTextContent);
                    end
                end
            end
        end
        PARAM.objFuncWeight.nbWeight = size(PARAM.objFuncWeight.value,2);
    end
           
end

