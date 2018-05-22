function randomizedlist = Randomizer(BoundariesFlowData,Method,flowside,InletMask,OutletMask,NumImagesPerCardiacCycle,NumberofHeartCycles,NumberofBoundaries,x,y,FolderNameList,InOutSwitch)

percenterror = input('What is the percent error to be applied to the sample data?: ');
error = rand(size(BoundariesFlowData))-0.5;
Simulated4DPCMRData = BoundariesFlowData + ((BoundariesFlowData).*error.*(percenterror/100));

InletList= zeros(NumImagesPerCardiacCycle,NumberofBoundaries);
OutletList=zeros(NumImagesPerCardiacCycle,NumberofBoundaries);
for i = 1:NumberofBoundaries
    if InOutSwitch(i) == 1
        InletList(:,i) = Simulated4DPCMRData(:,i);
    else
        OutletList(:,i) = Simulated4DPCMRData(:,i);
    end
end

sumOutlet = sum(OutletList,2);
sumInlet = sum(InletList,2);

switch Method
    case 1
        %% Flow rate correction(Ali's Method): Normalize and proportion the outlet or inlet flowrates to match the opposite side(whichever one is more trusted)
        % This method is good for when you have poor data on one side, most likely the outlet side
        
        switch flowside
            case 1
                BoundPercent = sumOutlet ./ sumInlet;
                NormalizeConstants = OutletMask.*(1./(BoundPercent));
            case 2
                BoundPercent = sumInlet ./ sumOutlet;
                NormalizeConstants = InletMask.*(1./(BoundPercent));
        end
        
        NormalizeConstants(NormalizeConstants == 0) = 1;
        BoundPercent(BoundPercent == 0) = 1;
        BoundAdjusted = Simulated4DPCMRData.*NormalizeConstants;
        
        %% Produce the final data for inlet/outlet comparison method
        
        BoundBCO1 = repmat(BoundAdjusted',1,NumberofHeartCycles);
        for i = 1:NumberofBoundaries
            randomizedlist(i,:) = interp1(y, BoundBCO1(i,:), x);
        end
    case 2
        %% All-boundaries flow-rate-proportional error correction (Alexa's Method)
        % This method is good for when you have good data and the outlet flow almost matches the inlet flow
        
        % Change the flowrate at each timestep such that the inlets match the outlets and the proportions remain the same
        error = sum(InletList,2)-sum(OutletList,2);
        
        % Find the subdivision of error with which to multiply each boundary flowrate
        subdiverror = error./(sum(Simulated4DPCMRData,2));
        
        % Generate a table of error corrections for each boundary condition at each timestep (negative on the inlet boundary to drop its flowrate)
        InletNegate = -InletMask;
        InletPiece = InletList.*InletNegate;
        proporterror = (InletPiece+OutletList).*subdiverror;
        
        % Add the error correction to the original values and generate a table of corrected flowrates that will be proportional to the original flowrate data
        FlowTable = proporterror+Simulated4DPCMRData;
        
        %% Produce the final data for proportional error method
        BoundBCO1 = repmat(FlowTable',1,NumberofHeartCycles);
        for i = 1:NumberofBoundaries
            randomizedlist(i,:) = interp1(y, BoundBCO1(i,:), x);
        end
end

figure;
hold on;
for i = 1:NumberofBoundaries
plot(x, randomizedlist(i,:));
end
xlabel('Time(s)');
ylabel('Mass flowrate (kg/s)');
legend([FolderNameList]);
title('Randomized Data')