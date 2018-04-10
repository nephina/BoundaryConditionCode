function randomizedlist = Randomizer(FlowTable,numofbounds,inoutswitch)

percent = input('What percentage should the error be limited to?: ');
randomMult = (rand(size(FlowTable)) - 0.5).*(percent/100);
RandomFlowTable = FlowTable + FlowTable.*randomMult;

inletList = zeros(size(FlowTable));
outletList = zeros(size(FlowTable));
for i = 1:numofbounds
    if inoutswitch(i) == 1
        inletList(:,i) = RandomFlowTable(:,i);
    else
        outletList(:,i) = RandomFlowTable(:,i);
    end
end

sumOutlet = sum(outletList,2);
sumInlet = sum(inletList,2);
OutletMask = outletList;
OutletMask(OutletMask > 0) = 1;
InletMask = inletList;
InletMask(InletMask > 0) = 1;

error = sum(inletList,2)-sum(outletList,2);
subdiverror = error./(sum(RandomFlowTable,2));
InletNegate = -InletMask;
InletPiece = inletList.*InletNegate;
proporterror = (InletPiece+outletList).*subdiverror;
randomizedlist = proporterror+RandomFlowTable;