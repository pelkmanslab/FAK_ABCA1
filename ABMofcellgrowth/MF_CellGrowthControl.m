function[CurrentModelMatrix, interactionmatrix] = MF_CellGrowthControl(energymatrix, CurrentModelMatrix, interactionmatrix,maxCellsize,minCellsize)
%[MF]Major modif on 06/13

%This module controls the growth of a cell.


%give a tag for border cells that should be frozen. Updated in the
%release_tension function.


if ~isfield(CurrentModelMatrix, 'FreezeTag')
    CurrentModelMatrix.FreezeTag = false(CurrentModelMatrix.numberofcells,1);
elseif size(CurrentModelMatrix.FreezeTag,1) < CurrentModelMatrix.numberofcells
    ExtraZeros = false(CurrentModelMatrix.numberofcells - size(CurrentModelMatrix.FreezeTag,1),1);
    CurrentModelMatrix.FreezeTag = [CurrentModelMatrix.FreezeTag;ExtraZeros];
end

for i = 1:CurrentModelMatrix.numberofcells
if CurrentModelMatrix.Nuclei_Location(i,1) <= 30 || CurrentModelMatrix.Nuclei_Location(i,2) <= 30 ...
        || CurrentModelMatrix.Nuclei_Location(i,1) >= CurrentModelMatrix.rownumber-30 ...
        || CurrentModelMatrix.Nuclei_Location(i,2) >= CurrentModelMatrix.columnnumber - 30
    CurrentModelMatrix.FreezeTag(i) = true;
end
end

%I take the surroundings of the cell to impact the nucleus growth
[nuclei_outlines,~]=MF_GetNucleiOutlinesFromVectData_img(CurrentModelMatrix);
CurrentModelMatrix.nuclei_area = MF_GetNucleiFromVectData_img(CurrentModelMatrix);
occupied_area=MF_GetCellArea_img(CurrentModelMatrix,interactionmatrix,true);
for iCell = 1:CurrentModelMatrix.numberofcells
energymatrix = MF_SetWholeCellPotential(CurrentModelMatrix, interactionmatrix, iCell,occupied_area, ...
    CurrentModelMatrix.nuclei_area,nuclei_outlines,energymatrix,1,true);%defines the accessible area for the cell
end

maxAccessibleArea = (CurrentModelMatrix.EnRadius*2+1)^2; %
minAccessibleArea =  ((sqrt(minCellsize./(pi))+5)^2)*pi; %The minimal area of a cell
accessible_surf = cellfun(@(x) sum(sum(x < 1)), energymatrix.StoredCellpotential, 'UniformOutput', true)';

RatioMaxMinSize = maxCellsize./minCellsize;%how much to multiply the min size to get the max size or to divide the max to get the min :) obviously, I am tired.
MaxDiff = log(maxAccessibleArea) - log(minAccessibleArea);
ActualDiff = log(maxAccessibleArea) - log(accessible_surf);%close to 0 if big, close to MaxDiff if small
%now let's modify the Ratio
ModifRatio = (RatioMaxMinSize-1).*(ActualDiff./MaxDiff);% = 1 when the actual diff is max, the modif ratio diminishes when the area get close to its max.
CellsizeLimit = round(maxCellsize./(1+ModifRatio));
GrowthVelocity = round(180./(1+ModifRatio));

%here is the rule to make cells growing
% t1(:,1)= exp(-CurrentModelMatrix.LCD)*CurrentModelMatrix.growthtimestepincrement.*6000;%
% t1(t1(:,1)>180.*CurrentModelMatrix.growthtimestepincrement,1)=180.*CurrentModelMatrix.growthtimestepincrement;
t1(:,1)=zeros(size(CurrentModelMatrix.LCD,1),1);t1(:,1) = GrowthVelocity;
ContinueGrowth = CurrentModelMatrix.Nuclei_Area < CellsizeLimit;

%Post mitotic and interphase cells grow in different ways.
CurrentModelMatrix.Nuclei_Area(CurrentModelMatrix.mitotictag==0 & ~CurrentModelMatrix.FreezeTag & ContinueGrowth)= ...
    CurrentModelMatrix.Nuclei_Area(CurrentModelMatrix.mitotictag==0 & ~CurrentModelMatrix.FreezeTag & ContinueGrowth) ...
    +t1(CurrentModelMatrix.mitotictag==0 & ~CurrentModelMatrix.FreezeTag & ContinueGrowth,1 );
CurrentModelMatrix.Nuclei_Area(CurrentModelMatrix.mitotictag~=0)=CurrentModelMatrix.Nuclei_Area(CurrentModelMatrix.mitotictag~=0)+1000*(log(CurrentModelMatrix.status(CurrentModelMatrix.mitotictag~=0)+1)-log(CurrentModelMatrix.status(CurrentModelMatrix.mitotictag~=0)));

%as nuclei are elipses, their definition is a little bit more complex than
%simple circle
occupied_area_PRIME=MF_GetCellArea_img(CurrentModelMatrix,interactionmatrix,true);
for iCell=1:CurrentModelMatrix.numberofcells
f = CurrentModelMatrix.radius(iCell,3);
% radius = zeros(CurrentModelMatrix.numberofcells,3);
oldradius=sqrt(CurrentModelMatrix.Nuclei_Area(iCell)./pi./(f.*(2-f)));
radius(1)=round(f.*oldradius);%new radius
radius(2)=round((1 +(1-f)).*oldradius);
radius(3)= f;
radius = sort(radius,2,'descend');

% bckup1 = CurrentModelMatrix.radius(:,1) ;
% bckup2 = CurrentModelMatrix.radius(:,2) ;
% bckup3 = CurrentModelMatrix.radius(:,3) ;

fprintf('Checking if growth of Nuclei %d is possible',iCell)
occupied_area = occupied_area_PRIME;
occupied_area(occupied_area == CurrentModelMatrix.ObjectID(iCell)) = 0;
occupied_area = logical(occupied_area);
CurrentModelMatrix.nuclei_area = MF_GetNucleiFromVectData_img(CurrentModelMatrix);
NucleiImage = CurrentModelMatrix.nuclei_area;
Image1 =  CurrentModelMatrix.ObjectID(iCell).* (NucleiImage == CurrentModelMatrix.ObjectID(iCell));
%imshow(Test1)
Test1 = (Image1 > 0) + occupied_area;
if any(Test1(:) == 2)
    disp(',...No!')
    CanItGrow(iCell) = false;
    continue
else
    CanItGrow(iCell) = true;
    disp(',...yes!')
CurrentModelMatrix.radius(iCell,1) = radius(1);
CurrentModelMatrix.radius(iCell,2) = radius(2);
CurrentModelMatrix.radius(iCell,3) = radius(3);
end

end
% %stop nuclei growth if growing on another cell cytoplasm

% for iCell=1:CurrentModelMatrix.numberofcells
% 
% 
% 
% % imshow(Image1)
% 
% count = 0;
% while any(Test1(:) == 2)
%     tic
%     count = count+1;
%     if count == 1
%         disp('a nuclei is over a close cell cytoplasm, stop its growth')
%     end 
%     
%     CurrentModelMatrix.radius(iCell,1) = bckup1(iCell);
%     CurrentModelMatrix.radius(iCell,2) = bckup2(iCell);
%     CurrentModelMatrix.radius(iCell,3) = bckup3(iCell);
%     
%     CurrentModelMatrix.nuclei_area = Display_nuclei_ellipsewithrotation_nodisplay(CurrentModelMatrix);
%     NucleiImage = CurrentModelMatrix.nuclei_area;
%     Image1 =  CurrentModelMatrix.ObjectID(iCell).* (NucleiImage == CurrentModelMatrix.ObjectID(iCell));
%     Test1 = (Image1 > 0) + occupied_area ;
%     
%     toc
% end
% 
% end


%get the list of overlapping nuclei

overlapping_list= MF_FindNucleiTension_list(CurrentModelMatrix);%, virtualimage2

%proceed with tension release as long as there are still overlapping nuclei
while ~isempty(overlapping_list) 
%release
CurrentModelMatrix=MF_ReleaseNucleiTensions(CurrentModelMatrix, overlapping_list);
%determine which nuclei overlap
overlapping_list= MF_FindNucleiTension_list(CurrentModelMatrix);%, virtualimage2)
%virtualimage2=Display_nuclei_ellipsewithrotation_nodisplay(ModelMatrix);
end




 

for iCell=1:CurrentModelMatrix.numberofcells
    if CurrentModelMatrix.FreezeTag(iCell) == false && CanItGrow(iCell) == true;
angle=interactionmatrix.bordernodes{iCell,1}(:,2);%all the node angles
%calculate here the distance from nuclei membrane to centroid  
factorl=zeros(length(angle),1);
nodeXpos=zeros(length(angle),1);
nodeYpos=zeros(length(angle),1);
for newnode=1:length(angle)
     
    factort=atan(CurrentModelMatrix.radius(iCell,1)/CurrentModelMatrix.radius(iCell,2)*tan(angle(newnode)-CurrentModelMatrix.angle(iCell,1)));
    factorl(newnode,1)=sqrt((CurrentModelMatrix.radius(iCell,1)*cos(factort))^2 + (CurrentModelMatrix.radius(iCell,2)*sin(factort))^2);%Norm calculation
    factorl(newnode,1)=factorl(newnode,1)+5;  
    dist1=factorl(newnode).*cos(angle(newnode));%Xshift
    dist2=factorl(newnode).*sin(angle(newnode));%Yshift

%now define its position in the image
nodeXpos(newnode,1)=round(CurrentModelMatrix.Nuclei_Location(iCell,1)+dist1);
nodeYpos(newnode,1)=round(CurrentModelMatrix.Nuclei_Location(iCell,2)+dist2);
end   
%is the distance enough? if yes update the data, if not go further in order
%to avoid node placements generating a border crossing a nuclei
index_correct=factorl>interactionmatrix.bordernodes{iCell,1}(:,1);
interactionmatrix.bordernodes{iCell,1}(index_correct,1)=factorl(index_correct,1);
interactionmatrix.bordernodes{iCell,1}(index_correct,3)=nodeXpos(index_correct,1);
interactionmatrix.bordernodes{iCell,1}(index_correct,4)=nodeYpos(index_correct,1);
interactionmatrix.bordernodes{iCell,1}(index_correct,6)=factorl(index_correct,1);

[~,interactionmatrix] = MF_BendMembrane(interactionmatrix,CurrentModelMatrix, iCell);
    end
 end


