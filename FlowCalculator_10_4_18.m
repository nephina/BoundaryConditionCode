% Before running, make sure to name all the boundary folders ending in "Boundary", and all the planes in the format
% plane1.0.0.csv,plane2.0.0.csv,etc. (in paraview you would write them to plane1.csv and it adds the ".0.0")
% Then check that all your timing variables are correct, and if you had the patience to get more than 3 planes, change that variable as well.

clc
clear all
close all

%% Set-up variables

global TimeStepOfFinalData;
global BeatsPerMin;
global NumImagesPerCardiacCycle;
global NumberofPlanes;
global NumberofHeartCycles;


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
        % BoundaryRawData(:,i,j) = abs(Bounds(1:end,7)); % This is a 3D matrix set up as (boundaryflow, sampleplane#, boundary#)
        BoundaryRawData(:,i,j) = Bounds(1:end,7); % No absolute value calculation, some boundaries might cross 0 flow
    end
end

%% Make all boundary datasets majorly positive (avoids the error of abs() for boundaries that may cross 0 flow)

for i = 1:NumberofBoundaries
    if mean(mean(BoundaryRawData(:,:,i))) >= 0
    else
        BoundaryRawData(:,:,i) = -BoundaryRawData(:,:,i);
    end
end

%% Flow conversion to kg/s and averaging to find the mean flowrate

for i = 1:NumberofBoundaries
    ConvertedFlowData = BoundaryRawData.*1060.*1e-9; % 1060*1e-9 is the conversion to mass for one cubic millimeter of blood
    BoundariesFlowData(:,i) = mean(ConvertedFlowData(:,:,i),2); % It then uses the mean function to combine the planes, so we have a 2D matrix of form (time, boundary)
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
InletMask = repmat(InOutSwitch',1,NumImagesPerCardiacCycle)';
OutletMask = 1 - InletMask;
FlipMask = ((-(1-InOutSwitch))); 
FlipMask(FlipMask == 0) = 1;

global deltaT;
global x;
global y;

deltaT = (60/BeatsPerMin)/NumImagesPerCardiacCycle;
y0 = 0:deltaT:((BeatsPerMin/60)-2*deltaT);
y = 0:deltaT:((7*((60/BeatsPerMin)))-(1/BeatsPerMin));
x = 0:TimeStepOfFinalData:y(end);

%% Plot the raw data and the averaged data of each boundary in order to make an informed choice on whether the data is good

RawDataDisplayTable = cat(3,BoundariesFlowData,permute(ConvertedFlowData,[1,3,2]));
RawDataDisplayTable = RawDataDisplayTable.*FlipMask;
plotbounds = [0 NumImagesPerCardiacCycle min(RawDataDisplayTable(:)) max(RawDataDisplayTable(:))]; % Scale each plot the same to get an accurate picture of how they compare
numofsubplots = (ceil(NumberofBoundaries/3)*3)+3; % Automatically figures out how many subplots are needed in a 2*n form to display all the boundaries without wasting space
subplotsydim = numofsubplots/3;

f = figure;
for i = 1:NumberofBoundaries
    hold on
    Color(i,:) = [((1/NumberofBoundaries)*i) (-((1/NumberofBoundaries)*i)+1) rand(1)]; % Trying to get a good differentiable spread of colors
    subplot(subplotsydim,3,i);
    % Display original .csv data
    for j = 1:NumberofPlanes
        hold on;
        plot(RawDataDisplayTable(:,i,(j+1)),'Color',Color(i,:),'--')
        axis(plotbounds);
        xlabel('Time in Heart Cycle');
        ylabel('Mass flowrate (kg/s)');
    end
    % Display mean of each boundary's .csv files
    plot(RawDataDisplayTable(:,i,1),'Color',Color(i,:),'Linewidth',1)
    title(BoundaryNameList(i));
    % Display means of .csv files as one cohesive plot
    subplot(subplotsydim,3,[numofsubplots-2 numofsubplots-1 numofsubplots]);
    plot(RawDataDisplayTable(:,i,1),'Color',Color(i,:),'Linewidth',1)
    axis(plotbounds);
    title('Mean Boundaries')
    xlabel('Time in Heart Cycle');
    ylabel('Mass flowrate (kg/s)');
end
% Function to let user adjust the plot before having to choose a method, I don't know why the pause(1) is needed but it makes it work in octave...
pause(1);
uiwait(gcf);

%% Adjustment methods for constant volume, equal inflow/outflow condition

Method = menu('What method to use?','Proportionized Using Reliable Flow','Proportionally Assigned Difference','quit');
close all;
flowside = 0;
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
        FlowTable = BoundariesFlowData.*NormalizeConstants;
        
        %% Produce the final data for inlet/outlet comparison method
        
        BoundBCO1 = repmat(FlowTable',1,NumberofHeartCycles);
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

%% Plotting

% Plot Output Data(all cycles)
for i = 1:NumberofBoundaries
    subplot(1,2,1);
    hold on;
    plot(x, BoundBCI(i,:),'Color',Color(i,:),'Linewidth',1);
end
% Plot single cycle, original data and recalculated data
for i = 1:NumberofBoundaries
    subplot(1,2,2);
    hold on;
    plot(FlipMask(i).*FlowTable(:,i),'Color',Color(i,:),'Linewidth',1);
    plot(RawDataDisplayTable(:,i,1),'Color',Color(i,:),'--');
end
hold on;
subplot(1,2,1);
xlabel('Time(s)');
ylabel('Mass flowrate (kg/s)');
title('Output Data');
subplot(1,2,2);
xlabel('Time in Heart Cycle');
ylabel('Mass flowrate (kg/s)');
title('Comparison To Original Data');

pause(1);
uiwait(gcf);

%% Adding randomness for Moji's code (takes the ground truth to be the generated boundary conditions, and overlays a simulated 4DPCMR error)

Randomness = menu('Add randomness to the boundary conditions?','no','yes');
if Randomness == 2
    BoundBCI = Randomizer(BoundariesFlowData,Method,flowside,InletMask,OutletMask,NumImagesPerCardiacCycle,NumberofHeartCycles,NumberofBoundaries,x,y,FolderNameList,InOutSwitch);
end

%% Output

app = menu('What are the boundary conditions for?','OpenFOAM','Ansys Fluent','No Output');
switch app
    case 1
       for i = 1:NumberofBoundaries
                    
           fileID = fopen(BoundaryNameList{i},'w');
           fprintf(fileID, '\t(\n');
           for j=1:length(x)
               fprintf(fileID, '\t(%d %d) \n', x(j), BoundBCI(i,j));
           end
           fprintf(fileID, '\t);\n');
           fclose(fileID);
        end
    case 2
        for i = 1:NumberofBoundaries
            fileName = strcat(BoundaryNameList{i},'.txt');
            fileID = fopen(fileName,'w');
            fprintf(fileID, ['((\n' BoundaryNameList{i} 'transient %d 1' '\n (time \n'], length(x));
            for j=1:length(x)
                fprintf(fileID, '\t %d \n', x(j));
            end
            fprintf(fileID, ' )\n (flowRate \n');
            for j=1:length(x)
                fprintf(fileID, '\t %d \n', BoundBCI(i,j));
            end
            fprintf(fileID, ' )\n)\n');
            fclose(fileID);
        end
    case 3
end