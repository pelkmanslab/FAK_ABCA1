function [virtualimage,CurrentModelMatrix,interactionmatrix,energymatrix] = MF_InitABmodel(StartNumOfCells,interval,totalenergy,NoN,radius)

load('C:\Users\MAT\Desktop\CurrentModelMat_init.mat');
load('C:\Users\MAT\Desktop\interactionmat_init.mat');
CurrentModelMatrix.growthtimestepincrement = interval;%this has to be said in Hours.
CurrentModelMatrix.rownumber = 3000;
CurrentModelMatrix.columnnumber = 3000;
interactionmatrix.numberofnodes = NoN;
%sample down the number of cells
NumOfCells = size(CurrentModelMatrix.Nuclei_Location,1);
Incr = floor(NumOfCells./StartNumOfCells);
PosTMP = 1:Incr:188;Pos2RMV = (1:188)';Pos2RMV = logical(Pos2RMV);Pos2RMV(PosTMP) = 0;
%go through the structures
FieldNamesCMM = fieldnames(CurrentModelMatrix);
FieldNamesIM  = fieldnames(interactionmatrix);
for i = 1:size(FieldNamesCMM,1)
    if size(CurrentModelMatrix.(FieldNamesCMM{i}),1) == NumOfCells && size(CurrentModelMatrix.(FieldNamesCMM{i}),2) < 5 %if it is a field to scale down
        CurrentModelMatrix.(FieldNamesCMM{i})(Pos2RMV,:) = []; 
    end
    if CurrentModelMatrix.(FieldNamesCMM{i})(1,1) == NumOfCells && size(CurrentModelMatrix.(FieldNamesCMM{i}),1) == 1 %if it is a field to scale down
        CurrentModelMatrix.(FieldNamesCMM{i}) = size(PosTMP,2); 
    end
end


for i = 1:size(FieldNamesIM,1)
    if size(interactionmatrix.(FieldNamesIM{i}),1) == NumOfCells && size(interactionmatrix.(FieldNamesIM{i}),2) < 5 %if it is a field to scale down
        interactionmatrix.(FieldNamesIM{i})(Pos2RMV,:) = []; 
    end
    if size(interactionmatrix.(FieldNamesIM{i}),1) == NumOfCells && size(interactionmatrix.(FieldNamesIM{i}),2) == NumOfCells %if it is a field to scale down
        interactionmatrix.(FieldNamesIM{i})(Pos2RMV,:) = [];
        interactionmatrix.(FieldNamesIM{i})(:,Pos2RMV) = [];
    end
end

Pos2RMV = CurrentModelMatrix.Nuclei_Location(:,1) <= 400 | CurrentModelMatrix.Nuclei_Location(:,2) <= 400;

for i = 1:size(FieldNamesCMM,1)
    if size(CurrentModelMatrix.(FieldNamesCMM{i}),1) > 1 && size(CurrentModelMatrix.(FieldNamesCMM{i}),2) < 5 %if it is a field to scale down
        CurrentModelMatrix.(FieldNamesCMM{i})(Pos2RMV,:) = []; 
    end
end
CurrentModelMatrix.numberofcells = size(Pos2RMV,1) - sum(Pos2RMV);
for i = 1:size(FieldNamesIM,1)
    if size(interactionmatrix.(FieldNamesIM{i}),1) > 1 && size(interactionmatrix.(FieldNamesIM{i}),2) < 5 %if it is a field to scale down
        interactionmatrix.(FieldNamesIM{i})(Pos2RMV,:) = []; 
    end
    if size(interactionmatrix.(FieldNamesIM{i}),1) == size(Pos2RMV,1) && size(interactionmatrix.(FieldNamesIM{i}),2) == size(Pos2RMV,1) %if it is a field to scale down
        interactionmatrix.(FieldNamesIM{i})(Pos2RMV,:) = [];
        interactionmatrix.(FieldNamesIM{i})(:,Pos2RMV) = [];
    end
end


for i = 1:CurrentModelMatrix.numberofcells
    [interactionmatrix] = MF_ResetNodePosition(interactionmatrix,CurrentModelMatrix, i,20);
end
% for i = 1:size(FieldNamesCMM,1) 
%     
% end
% 
% for i = 1:size(FieldNamesIM,1)
% end

%randomly attributes a point in the cell cycle.
TOG4CellCycle = 22 .* rand(CurrentModelMatrix.numberofcells,1); %we assume from movies that a normal cell cycle is 22 hours, variability comes afterwards.
CurrentModelMatrix.status(:) = TOG4CellCycle;
CurrentModelMatrix.AbsolutePos = 1:CurrentModelMatrix.columnnumber*CurrentModelMatrix.rownumber;
CurrentModelMatrix.AbsolutePos = uint32(reshape(CurrentModelMatrix.AbsolutePos,CurrentModelMatrix.rownumber,CurrentModelMatrix.columnnumber));
virtualimage = MF_ModelBckbone_Img(interactionmatrix,CurrentModelMatrix);

%initialization of parameters
energymatrix.factors = [1,1,1];
energymatrix.factors_feature={'multiplication factor for whole cell potential','factor for local (distance) potential','factor for conservativness of cell movement'};
energymatrix.totalenergy=zeros(CurrentModelMatrix.numberofcells,1);
%minimal_dist_to_nucleimembrane=3;

%energymatrix initialization
for i=1:CurrentModelMatrix.numberofcells
    energymatrix.nodeenergy{i,1}=zeros(interactionmatrix.numberofnodes,1);
    energymatrix.distributionOfenergy{i,1}=zeros(interactionmatrix.numberofnodes,1);
    energymatrix.totalenergy(i,1)=random('normal',totalenergy(1,1),totalenergy(1,2));%
    energymatrix.distributionOfenergy{i}(:,1) = energymatrix.totalenergy(i,1)/NoN; %init, for the moment equal distrib of energy
    energymatrix.cell_vector(i,:)=[0,0];
end

for i=1:CurrentModelMatrix.numberofcells
interactionmatrix.bordernodes{i,1}(:,7)=zeros(interactionmatrix.numberofbordernodes(i),1);
end

% factor = 4;
% ratio = factor/interactionmatrix.numberofnodes;
% DistribOfE = ((1/factor):ratio:factor)';%linear distribution, let see.
% 
% if size(DistribOfE,1)>interactionmatrix.numberofnodes
% DistribOfE=DistribOfE(1:interactionmatrix.numberofnodes);
% elseif size(DistribOfE,1)<interactionmatrix.numberofnodes
% newPos = interactionmatrix.numberofnodes - size(DistribOfE,1);
% newPos = (1:newPos)'; 
% newVal(newPos) = DistribOfE(end) + (newPos .* ratio);
% newPos = newPos + size(DistribOfE,1);
% DistribOfE(newPos,1) = newVal;
% end
% interactionmatrix.DistribOfE = flipud(DistribOfE);
DistribOfE = ones(interactionmatrix.numberofnodes,1);

 interactionmatrix.DistribOfE = DistribOfE;

%lets draw a generic energy landscape:

CurrentModelMatrix.EnRadius = radius;
% a = (0.01:1/100:1); a = [a sort(a,'descend')];
% b = a';
% for i = 1:200;
% aArr(i,:) = a;
% bArr(:,i) = b;
% end
% Ldscp = aArr.*bArr;

Image = zeros(radius.*2+1); Image(radius+1,radius+1) = 1;
ImageFinal = Image;
for i = 1:2:(radius-round(radius/4))

SE = strel('disk',i+1);
Image2 = imdilate(Image,SE);
ImageFinal = ImageFinal + Image2;

end
ImageFinal = ImageFinal./max(ImageFinal(:));
CurrentModelMatrix.basicLandscape = ImageFinal;


% CurrentModelMatrix.smoothNucleiEnergy = ones(251,251);
% CurrentModelMatrix.smoothNucleiEnergy = Gauss2D(CurrentModelMatrix.smoothNucleiEnergy,50);

end

