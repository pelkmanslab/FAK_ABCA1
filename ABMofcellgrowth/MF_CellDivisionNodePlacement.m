function[CurrentModelMatrix,energymatrix,virtualimage]=MF_CellDivisionNodePlacement(CurrentModelMatrix, virtualimage, allcells, minCellsize,interactionmatrix,energymatrix)
%[MF] massive changes 06/13
% CELLDIVISION(CurrentModelMatrix,virtualimage, whichcelldivides)
% Makes out of old nucleus,two new nuclei. Hence it simulates the process
% of cell division and updates CurrentModelMatrix (but not Nucleus_Area!
% and LCD measurement!

%for redo
Rerun1 =CurrentModelMatrix;
Rerun2 =virtualimage;
Rerun3 =allcells;
Rerun4 =minCellsize;
Rerun5 =interactionmatrix;
mkspace = 0;
mark = true;
while mark == true
    mark = false;
    
CurrentModelMatrix = Rerun1; 
virtualimage = Rerun2;
allcells = Rerun3;
minCellsize = Rerun4;
interactionmatrix = Rerun5;  



occupied_area_PRIME=MF_GetCellArea_img(CurrentModelMatrix,interactionmatrix,true);


if ~isempty(allcells)
for whichcelldivides = allcells %there might be more than one cell dividing per frame (loop over all those)
 if mkspace > 0
     trshld = 30+(mkspace.*4);
    [CurrentModelMatrix,interactionmatrix,energymatrix] = MF_MakeSpace4DividingCells(CurrentModelMatrix,interactionmatrix,energymatrix,whichcelldivides,trshld); 
 end
    
oldID = CurrentModelMatrix.ObjectID(whichcelldivides);
%First: How many cells are there? => where to place new cells?
numberofcells = CurrentModelMatrix.numberofcells;
%Location of new cells, also save new angle
%ne could think about coding here the fact that division axis might be
%constrained by the islet of cells
oldcoordinates = CurrentModelMatrix.Nuclei_Location(whichcelldivides,:);
phi = 2 * pi *random('Uniform',0,1);
OffSetV = random('normal',20,3);
CurrentModelMatrix.Nuclei_Location(whichcelldivides,:) = oldcoordinates + (CurrentModelMatrix.radius(whichcelldivides)+OffSetV)/2* [cos(phi) sin(phi)];
CurrentModelMatrix.Nuclei_Location(numberofcells+1,:) = oldcoordinates - (CurrentModelMatrix.radius(whichcelldivides)+OffSetV)/2* [cos(phi) sin(phi)];

%let's try this, it should be quite efficient for fixing border cells cases
if CurrentModelMatrix.Nuclei_Location(whichcelldivides,1) < 10 || ...
    CurrentModelMatrix.Nuclei_Location(whichcelldivides,1) > CurrentModelMatrix.rownumber-10 || ...
    CurrentModelMatrix.Nuclei_Location(whichcelldivides,2) < 10 || ...
    CurrentModelMatrix.Nuclei_Location(whichcelldivides,2) > CurrentModelMatrix.columnnumber-10
CurrentModelMatrix.FreezeTag(whichcelldivides) = true;
end
if CurrentModelMatrix.Nuclei_Location(numberofcells+1,1) < 10 || ...
    CurrentModelMatrix.Nuclei_Location(numberofcells+1,1) > CurrentModelMatrix.rownumber-10 || ...
    CurrentModelMatrix.Nuclei_Location(numberofcells+1,2) < 10 || ...
    CurrentModelMatrix.Nuclei_Location(numberofcells+1,2) > CurrentModelMatrix.columnnumber-10
CurrentModelMatrix.FreezeTag(numberofcells+1) = true;
else
CurrentModelMatrix.FreezeTag(numberofcells+1) = false;    
end
%


CurrentModelMatrix.angle(whichcelldivides,:) = phi + pi/2;
CurrentModelMatrix.angle(numberofcells+1,:) = phi + pi/2;

%Radius of new cells
% radius=zeros(2,3);
 fact=random('uniform',0.9,1.1);
%  f = CurrentModelMatrix.radius(whichcelldivides,3);
% Ratio = random('normal',0.5,0.05);
% oldradius = sqrt(CurrentModelMatrix.Nuclei_Area(whichcelldivides)/Ratio/pi/(f*(2-f)));%sqrt(CurrentModelMatrix.Nuclei_Area(whichcelldivides)/3/pi);
% radius(:,1)=round(fact.*oldradius);
% radius(:,2)=round((1 +(1-fact)).*oldradius(:,1));
% radius(:,3)= fact;
% radius = sort(radius,2,'descend');
% CurrentModelMatrix.radius(whichcelldivides,:) = radius(1,:);
% CurrentModelMatrix.radius(numberofcells+1,:) = radius(2,:);
% %Area of new cells
% CurrentModelMatrix.Nuclei_Area(whichcelldivides)= pi * radius(1,1) * radius(1,2);
% CurrentModelMatrix.Nuclei_Area(numberofcells+1) = pi * radius(2,1) * radius(2,2);
% %set minimal area of 1000
% if CurrentModelMatrix.Nuclei_Area(whichcelldivides) < 1000 || ...
%         CurrentModelMatrix.Nuclei_Area(numberofcells+1) < 1000
CurrentModelMatrix.Nuclei_Area(whichcelldivides) = minCellsize;
CurrentModelMatrix.Nuclei_Area(numberofcells+1) = minCellsize;

oldradius = sqrt(CurrentModelMatrix.Nuclei_Area(whichcelldivides)/pi/(fact*(2-fact)));
radius(:,1) = round(fact.*oldradius);
radius(:,2)= round((1 +(1-fact)).*oldradius(:,1));
radius(:,3)= fact;
CurrentModelMatrix.radius(whichcelldivides,:) = radius(1,:);

oldradius = sqrt(CurrentModelMatrix.Nuclei_Area(numberofcells+1)/pi/(fact*(2-fact)));
radius(:,1) = round(fact.*oldradius);
radius(:,2)= round((1 +(1-fact)).*oldradius(:,1));
radius(:,3)= fact;
CurrentModelMatrix.radius(numberofcells+1,:) = radius(1,:);



%delete old nuclei in image
virtualimage(virtualimage == CurrentModelMatrix.ObjectID(whichcelldivides)) = 0;
%Update of ObjectID, ParentID, FamilyID
% if max(CurrentModelMatrixD.ObjectID) + 1 ==400
%     newlabel=420;
% else
    newlabel= max(CurrentModelMatrix.ObjectID) + 1;
% end
CurrentModelMatrix.ObjectID(whichcelldivides)=newlabel;
CurrentModelMatrix.ObjectID(numberofcells+1)=newlabel+1;
CurrentModelMatrix.ParentID(numberofcells+1)=CurrentModelMatrix.ObjectID(whichcelldivides);
CurrentModelMatrix.FamilyID(numberofcells+1)=CurrentModelMatrix.FamilyID(whichcelldivides);
%Update of Color
CurrentModelMatrix.color(numberofcells+1,:)=CurrentModelMatrix.color(whichcelldivides,:);
%initialize new pFAK content of new cell
CurrentModelMatrix.pFAK4myObjects(numberofcells+1) = CurrentModelMatrix.pFAK4myObjects(whichcelldivides);
%Update Cell status
CurrentModelMatrix.status(numberofcells+1)=0 + random('Uniform',0,CurrentModelMatrix.growthtimestepincrement);%actually you don't know when occured the division between the two frames
CurrentModelMatrix.status(whichcelldivides)=0 + random('Uniform',0,CurrentModelMatrix.growthtimestepincrement);%same
%Update at last number of cells
CurrentModelMatrix.numberofcells = numberofcells + 1;
%set mitotic tag
CurrentModelMatrix.mitotictag(whichcelldivides) = 1;%why the hell do you put a 5, keep it logical
CurrentModelMatrix.mitotictag(numberofcells+1) = 1;

%safetycheck, do not put a nuclei on another cell node
%first create the good image
occupied_area = occupied_area_PRIME;
occupied_area(occupied_area == oldID) = 0;
occupied_area = logical(occupied_area);
CurrentModelMatrix.nuclei_area = MF_GetNucleiFromVectData_img(CurrentModelMatrix);
NucleiImage = CurrentModelMatrix.nuclei_area;
Image1 =  CurrentModelMatrix.ObjectID(whichcelldivides).* (NucleiImage == CurrentModelMatrix.ObjectID(whichcelldivides));
Image2 = CurrentModelMatrix.ObjectID(numberofcells+1).*(NucleiImage == CurrentModelMatrix.ObjectID(numberofcells+1));

% imshow(Image1)
Test1 = (Image1 > 0) + occupied_area;
Test2 = (Image2 > 0) + occupied_area;

oldLoc1 = CurrentModelMatrix.Nuclei_Location(whichcelldivides,:);
count = 0;
while any(Test1(:) == 2)
    tic
    count = count+1;
    if count == 1
        disp('a new nuclei is over a close cell cytoplasm, fixing')
    end
    CurrentModelMatrix.Nuclei_Location(whichcelldivides,:) = oldLoc1;
    
    if count < 10
        std = 1;
    else
        mark = true;
        mkspace = mkspace+1;
        disp('this division is not possible, rerun')        
        break
    end
    
    shiftX = round(random('normal',0,std));
    shiftY = round(random('normal',0,std));
    CurrentModelMatrix.Nuclei_Location(whichcelldivides,:) = CurrentModelMatrix.Nuclei_Location(whichcelldivides,:) + [shiftX shiftY];
    CurrentModelMatrix.nuclei_area = MF_GetNucleiFromVectData_img(CurrentModelMatrix);
    NucleiImage = CurrentModelMatrix.nuclei_area;
    Image1 =  CurrentModelMatrix.ObjectID(whichcelldivides).* (NucleiImage == CurrentModelMatrix.ObjectID(whichcelldivides));
    Test1 = (Image1 > 0) + occupied_area + (Image2 > 0);
    toc
end

oldLoc2 = CurrentModelMatrix.Nuclei_Location(numberofcells+1,:);
count = 0;
while any(Test2(:) == 2) && mark == false
    tic
    count = count+1;
    if count == 1
        disp('a new nuclei is over a close cell cytoplasm, fixing')
    end
    CurrentModelMatrix.Nuclei_Location(numberofcells+1,:) = oldLoc2;
    
    if count < 10
        std = 1;
    else
        mark = true;
        mkspace = mkspace+1;
        disp('this division is not possible, rerun')
        break
    end
    
    shiftX = round(random('normal',0,std));
    shiftY = round(random('normal',0,std));
    CurrentModelMatrix.Nuclei_Location(numberofcells+1,:) = CurrentModelMatrix.Nuclei_Location(numberofcells+1,:) + [shiftX shiftY];
    CurrentModelMatrix.nuclei_area = MF_GetNucleiFromVectData_img(CurrentModelMatrix);
    NucleiImage = CurrentModelMatrix.nuclei_area;
    Image2 = CurrentModelMatrix.ObjectID(numberofcells+1).*(NucleiImage == CurrentModelMatrix.ObjectID(numberofcells+1));
    Test2 = (Image2 > 0) + occupied_area + (Image1 > 0);
    toc

end

end
end
end
