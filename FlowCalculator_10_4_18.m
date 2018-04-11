% Before running, make sure to name all the boundary folders ending in "Boundary", and all the planes in the format
% plane1.0.0.csv,plane2.0.0.csv,etc. (in paraview you would write them to plane1.csv and it adds the ".0.0")
% Then check that all your timing variables are correct, and if you had the patience to get more than 3 planes, change that variable as well.

clc
clear all
close all

%% Set-up variables

TimeStepOfFinalData = 0.01;
BeatsPerMin = 87;
NumImagesPerCardiacCycle = 20;
NumberofPlanes = 3;
NumberofHeartCycles = 7;

%--------------------------------------- BEGIN CODE ---------------------------------------------------------------------------------------

%% Search the directory for all boundary folders with "Boundary" in their names and set up naming for the boundary files

BoundarySearch = dir('*Boundary*');
NumberofBoundaries = size(BoundarySearch);
NumberofBoundaries = NumberofBoundaries(1);
BoundaryTable = struct2cell(BoundarySearch);
BoundaryNameList= {};
for i = 1:NumberofBoundaries
    BoundaryNames = BoundaryTable(1,i,1);
    BoundaryName = BoundaryNames{1};
    FolderNameList{i} = BoundaryName;
    BoundaryName = BoundaryName(1:end-8);
    BoundaryNameList{i} = BoundaryName;
end

%% Read .csv data from each plane for each boundary
for j = 1:NumberofBoundaries
    for i = 1:NumberofPlanes;
        planeName = sprintf('plane%d0.0.csv', i);
        loadBounds =  sprintf('%s/%s', FolderNameList{j}, planeName);
        Bounds = dlmread(loadBounds,",", 1, 0);
        BoundaryRawData(:,i,j) = abs(Bounds(1:end,7)); % This is a 3D matrix set up as (boundaryflow, sampleplane#, boundary#)
    end
end

%% Flow conversion to kg/s and averaging to find the mean flowrate

for i = 1:NumberofBoundaries
    BoundariesFlowData(:,i) = mean(BoundaryRawData(:,:,i),2).*1060* 1e-9; % 1060*1e-9 is the conversion to mass for one cubic millimeter of blood
        % It then uses the mean function to combine the planes, so we have a 2D matrix of form (boundary, time)
end

%% Setup of outlets/inlets using user input

for i = 1:NumberofBoundaries
    inletoutletquery = sprintf('Is %s an inlet? (use 1/0):  ',BoundaryNameList{i});
    InOutSwitch(i) = input(inletoutletquery);
    clc
end

InletList= zeros(NumImagesPerCardiacCycle,NumberofBoundaries);
OutletList=zeros(NumImagesPerCardiacCycle,NumberofBoundaries);
for i = 1:NumberofBoundaries
    if InOutSwitch(i) == 1
        InletList(:,i) = BoundariesFlowData(:,i);
    else
        OutletList(:,i) = BoundariesFlowData(:,i);
    end
end

sumOutlet = sum(OutletList,2);
sumInlet = sum(InletList,2);
OutletMask = OutletList;
OutletMask(OutletMask > 0) = 1;
InletMask = InletList;
InletMask(InletMask > 0) = 1;

deltaT = (60/BeatsPerMin)/NumImagesPerCardiacCycle;
y0 = 0:deltaT:((BeatsPerMin/60)-2*deltaT);
y = 0:deltaT:((7*((60/BeatsPerMin)))-(1/BeatsPerMin));
x = 0:TimeStepOfFinalData:y(end);

%% Display the current boundary data and select a correction method

figure;
hold on;
for i = 1:NumberofBoundaries
    plot(1:NumImagesPerCardiacCycle, BoundariesFlowData(:,i));
end
xlabel('Time(s)');
ylabel('Mass flowrate (kg/s)');
legend([FolderNameList]);

Method = menu('What method to use?','Outlets/Inlets Proportionized Based on Most Reliable Flow','Proportionally Assigned Difference','quit');
close all;
switch Method
    case 1
        %% Flow rate correction(Ali's Method): Normalize and proportion the outlet or inlet flowrates to match the opposite side(whichever one is more trusted)
        % This method is good for when you have poor data on one side, most likely the outlet side
        
        flowside = menu('Which is more reasonable?','Inlet Flow Rates','Outlet Flow Rates');
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
        BoundAdjusted = BoundariesFlowData.*NormalizeConstants;
        
        %% Produce the final data for inlet/outlet comparison method
        
        BoundBCO1 = repmat(BoundAdjusted',1,NumberofHeartCycles);
        for i = 1:NumberofBoundaries
            BoundBCI(i,:) = interp1(y, BoundBCO1(i,:), x);
        end
    case 2
        %% All-boundaries flow-rate-proportional error correction (Alexa's Method)
        % This method is good for when you have good data and the outlet flow almost matches the inlet flow
        
        % Change the flowrate at each timestep such that the inlets match the outlets and the proportions remain the same
        error = sum(InletList,2)-sum(OutletList,2);
        
        % Find the subdivision of error with which to multiply each boundary flowrate
        subdiverror = error./(sum(BoundariesFlowData,2));
        
        % Generate a table of error corrections for each boundary condition at each timestep (negative on the inlet boundary to drop its flowrate)
        InletNegate = -InletMask;
        InletPiece = InletList.*InletNegate;
        proporterror = (InletPiece+OutletList).*subdiverror;
        
        % Add the error correction to the original values and generate a table of corrected flowrates that will be proportional to the original flowrate data
        FlowTable = proporterror+BoundariesFlowData;
        
        %% Produce the final data for proportional error method
        BoundBCO1 = repmat(FlowTable',1,NumberofHeartCycles);
        for i = 1:NumberofBoundaries
            BoundBCI(i,:) = interp1(y, BoundBCO1(i,:), x);
        end
end

%% Adding randomness for Moji's code (takes the ground truth to be the generated boundary conditions, and overlays a simulated 4DPCMR error)

Randomness = menu('Add randomness to the boundary conditions?','yes','no');
if Randomness == 1
    BoundBCI = Randomizer(BoundBCI,InletMask,InOutSwitch);
end

%% Plotting
figure;
hold on;
for i = 1:NumberofBoundaries
plot(x, BoundBCI(i,:));
end
xlabel('Time(s)');
ylabel('Mass flowrate (kg/s)');
legend([FolderNameList]);

%% Output

app = menu('What are the boundary conditions for?','OpenFOAM','Ansys Fluent','quit');
switch app
    case 1
        switch Method
            case 1
                for i = 1:NumberofBoundaries
                    
                    fileID = fopen(BoundaryNameList{i},'w');
                    fprintf(fileID, '(\n');
                    for j=1:length(x)
                        fprintf(fileID, '\t(%d %d) \n', x(j), BoundBCI(i,j));
                    end
                    fprintf(fileID, ');\n');
                    fclose(fileID);
                end
            case 2
                
                for i = 1:NumberofBoundaries
                    
                    fileID = fopen(BoundaryNameList{i},'w');
                    fprintf(fileID, '(\n');
                    for j=1:length(x)
                        fprintf(fileID, '\t(%d %d) \n', x(j), BoundBCI(i,j));
                    end
                    fprintf(fileID, ');\n');
                    fclose(fileID);
                end
            case 3
        end
    case 2
        switch Method
            case 1
                for i = 1:NumberofBoundaries
                    
                    fileID = fopen(BoundaryNameList{i},'w');
                    fprintf(fileID, ['(\n' BoundaryNameList{i} '\n (time \n'], length(x));
                    for j=1:length(x)
                        fprintf(fileID, '\t %d \n', x(j));
                    end
                    fprintf(fileID, ' )\n (flowRate \n');
                    for j=1:length(x)
                        fprintf(fileID, '\t %d \n', BoundBCI(j));
                    end
                    fprintf(fileID, ' )\n)\n');
                    fclose(fileID);
                end
            case 2
                for i = 1:NumberofBoundaries
                    
                    fileID = fopen(BoundaryNameList{i},'w');
                    fprintf(fileID, ['(\n' BoundaryNameList{i} '\n (time \n'], length(x));
                    for j=1:length(x)
                        fprintf(fileID, '\t %d \n', x(j));
                    end
                    fprintf(fileID, ' )\n (flowRate \n');
                    for j=1:length(x)
                        fprintf(fileID, '\t %d \n', BoundBCI(j));
                    end
                    fprintf(fileID, ' )\n)\n');
                    fclose(fileID);
                end
            case 3
        end
    case 3
end
