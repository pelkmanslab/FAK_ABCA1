function[temp4, interactionmatrix, labeled_overlapregions]=MF_FindOverlapsAndGuiltyNodes(CurrentModelMatrix, interactionmatrix)%, energymatrix)
%[MF] massive changes, the whole part of checking pixels in or not can be
%changed, the logic is not optimized I think. 06/13, work on going



%here is how it works: first it checks for each cells, what nodes are fine,
%inside the image. The module draw then the outlines made of segments
%between nodes (stored counterclockwise), finds the regions of the image
%corresponding to the outlines or to the interrior of all cells and defines
%the overlaps between cells, identifying coordinates corresponding to the
%intruder cell


%let's store this for later
PosInsideCell = cell(CurrentModelMatrix.numberofcells,1);%positions of the current image that are in the object concerned
AbsolutePos = cell(CurrentModelMatrix.numberofcells,1);%the actual patch of the object of interrest, with absolute pos in the total image > fast indexing
IMGCellpatchesFull = cell(CurrentModelMatrix.numberofcells,1);
BoundariesClosed = cell(CurrentModelMatrix.numberofcells,1);
% IdOverlapPix4iCell = cell(CurrentModelMatrix.numberofcells,1);%same patch always, but with only overlapping pixels involving the current object of interrest
% IMGCellpatches = cell(CurrentModelMatrix.numberofcells,1);%the actual patch of the object of interrest, contains the object boundaries as objID pixels 
% BoundariestoFill = cell(CurrentModelMatrix.numberofcells,1);%same patch always, with outlines in a logical array

interactionmatrix.nodestochange=[];
interactionmatrix.reservenode=cell(CurrentModelMatrix.numberofcells,CurrentModelMatrix.numberofcells);
interactionmatrix.extremaList=cell(CurrentModelMatrix.numberofcells,CurrentModelMatrix.numberofcells);
for i=1:CurrentModelMatrix.numberofcells
interactionmatrix.nodestochange{i,1}=zeros(interactionmatrix.numberofbordernodes(i),1);
end
for i=1:CurrentModelMatrix.numberofcells
interactionmatrix.bordernodes{i,1}(:,7)=zeros(interactionmatrix.numberofbordernodes(i),1);
end

nodesofcell=interactionmatrix.bordernodes;


%make an image with all boundaries
overlapimage=zeros(CurrentModelMatrix.rownumber, CurrentModelMatrix.columnnumber);

%totalindices = zeros(560.*CurrentModelMatrix.numberofcells,1);%7 nodes * 80 as maxi line length
totalindices = cell(CurrentModelMatrix.numberofcells,1);
totalindicesCELLID = cell(CurrentModelMatrix.numberofcells,1);
box_coord=zeros(CurrentModelMatrix.numberofcells,4);
node_list = cell(CurrentModelMatrix.numberofcells,1);
CoorXobj = cell(CurrentModelMatrix.numberofcells,1);
CoorYobj = cell(CurrentModelMatrix.numberofcells,1);
CoorYRel = cell(CurrentModelMatrix.numberofcells,1);
CoorXRel = cell(CurrentModelMatrix.numberofcells,1);

for whichcell=1:CurrentModelMatrix.numberofcells
    node_list{whichcell,1}=false(interactionmatrix.numberofbordernodes(whichcell),1);
    if ~CurrentModelMatrix.FreezeTag(whichcell)
%[MF] check for the eligibility of the node, checks if in image in
%particular.
for inde=1:interactionmatrix.numberofbordernodes(whichcell)
   if (nodesofcell{whichcell}(inde,3)>=1 && nodesofcell{whichcell}(inde,4)>=1 ...
       && nodesofcell{whichcell}(inde,4)<= CurrentModelMatrix.rownumber ... 
       && nodesofcell{whichcell}(inde,3)<= CurrentModelMatrix.columnnumber)
   %if the node is in the image we can process it, if not freeze the cell.
   %we can come on this later
   node_list{whichcell,1}(inde,1)=true; %if the node can be treated, it gets a number
   else
       CurrentModelMatrix.FreezeTag(whichcell) = true;
       %at this point the other nodes will be treated, at the next round
       %the cell won't be looked at anymore.
   end
end
% node_list{whichcell,1}(node_list{whichcell,1} == 0) = [];
% %[MF]we have now the nodes that can be treated 
% node_list{whichcell,1}=node_list{whichcell,1}';
%  coordOfNodesInCH = [];
%  TMP = nodesofcell{whichcell}(node_list{whichcell,1},3:4);coordOfNodesInCH(:,1) = TMP(:,2);coordOfNodesInCH(:,2) = TMP(:,1);
%  SelectedNodes = convhull(coordOfNodesInCH); SelectedNodes(end) = []; %we will gain a lot of time, however nodes that are inside the hull are not shown
%  node_list{whichcell,1} = node_list{whichcell,1}(SelectedNodes);
%  clear SelectedNodes TMP

leftboundary=min(interactionmatrix.bordernodes{whichcell,1}(node_list{whichcell,1},3));
rightboundary=max(interactionmatrix.bordernodes{whichcell,1}(node_list{whichcell,1},3));
topboundary=min(interactionmatrix.bordernodes{whichcell,1}(node_list{whichcell,1},4));
lowerboundary=max(interactionmatrix.bordernodes{whichcell,1}(node_list{whichcell,1},4));
box_coord(whichcell,:)=[topboundary,lowerboundary,leftboundary,rightboundary];


% currentIX = 0;
%let draw lines between nodes...
indices = cell(size(node_list{whichcell,1},1),1);
Node2LookAtX = nodesofcell{whichcell}(node_list{whichcell,1},4);
Node2LookAtY = nodesofcell{whichcell}(node_list{whichcell,1},3);
ParForList = (1:size(Node2LookAtX,1))' ;
FirstVal = ParForList(1,1);
lastVal = ParForList(end,1);
% Node2LookAtXRel{whichcell} = Node2LookAtX - (box_coord(whichcell,1)-1);%minus 1 of course as a coordinate equal to 0 is a nonsens
% Node2LookAtYRel{whichcell} = Node2LookAtY - (box_coord(whichcell,3)-1);
% Node2LookAtXRel{whichcell}(Node2LookAtXRel{whichcell} <= 0) = 1;
% Node2LookAtYRel{whichcell}(Node2LookAtYRel{whichcell} <= 0) = 1;
for i=1:size(Node2LookAtX,1)
   inde = ParForList(i,1);   
    if inde==lastVal %for the last node in the list
          newindices=MF_DrawLinesBetwnNodes_coord(Node2LookAtX(inde), ... 
          Node2LookAtY(inde),Node2LookAtX(FirstVal),Node2LookAtY(FirstVal));
          indices{i} = newindices ;
%         bndrTMP = size(newindices,1);
%         indices(currentIX+1:currentIX+bndrTMP,1) = newindices ;
%         currentIX = bndrTMP + currentIX;            
    else %for the others
          indePlus1 = ParForList(i+1,1);
          newindices=MF_DrawLinesBetwnNodes_coord(Node2LookAtX(inde), ... 
          Node2LookAtY(inde),Node2LookAtX(indePlus1),Node2LookAtY(indePlus1));
          indices{i} = newindices ;
%         bndrTMP = size(newindices,1);        
%         indices(currentIX+1:currentIX+bndrTMP,1) = newindices ; %indices of points to reconstitute boundaries of cells
%         currentIX = bndrTMP + currentIX;
    end
end
CoorXobj{whichcell} = cat(1,indices{:});%there can be redundant indices, but this is no problem
CoorYobj{whichcell} = CoorXobj{whichcell}(:,2);
CoorXobj{whichcell}(:,2) = [];
clear indices
%TMPimg = false(max(Node2LookAtXRel),max(Node2LookAtYRel));
% ShowPixels{whichcell} = false(max(Node2LookAtXRel{whichcell})+20,max(Node2LookAtYRel{whichcell})+20);
% for iii = 1:size(Node2LookAtXRel{whichcell},1)
%     ShowPixels{whichcell}(Node2LookAtXRel{whichcell}(iii),Node2LookAtYRel{whichcell}(iii)) = true;
% end
% IMGCellpatchesFull{whichcell} = poly2mask(Node2LookAtXRel{whichcell},Node2LookAtYRel{whichcell},max(Node2LookAtXRel{whichcell})+20,max(Node2LookAtYRel{whichcell})+20);
AbsolutePos{whichcell} = CurrentModelMatrix.AbsolutePos(box_coord(whichcell,1):box_coord(whichcell,2),box_coord(whichcell,3):box_coord(whichcell,4));
CoorXRel{whichcell} = CoorXobj{whichcell} - (box_coord(whichcell,1)-1);%minus 1 of course as a coordinate equal to 0 is a nonsens
CoorYRel{whichcell} = CoorYobj{whichcell} - (box_coord(whichcell,3)-1);
CoorXRel{whichcell}(CoorXRel{whichcell} <= 0) = 1;
CoorYRel{whichcell}(CoorYRel{whichcell} <= 0) = 1;
relPos = ((CoorYRel{whichcell}-1).*size(AbsolutePos{whichcell},1))+ CoorXRel{whichcell};
IMGCellpatchesFull{whichcell} = false(size(AbsolutePos{whichcell}));
IMGCellpatchesFull{whichcell}(relPos) = true;

BoundariesClosed{whichcell} = bwmorph(IMGCellpatchesFull{whichcell},'bridge');
% BoundariesClosed{whichcell} = bwmorph(BoundariesClosed{whichcell},'remove');
IMGCellpatchesFull{whichcell} = imfill(BoundariesClosed{whichcell},'holes');

indices = AbsolutePos{whichcell}(BoundariesClosed{whichcell});
PosInsideCell{whichcell} = AbsolutePos{whichcell}(IMGCellpatchesFull{whichcell});
indices(indices == 0) = []; %actually this should not occur.
% indices = indices([true;diff(indices(:))>0]);%keep in mind this replaces unique
% indices = indices([true;diff(indices(:))>0]);%because there can be overlapping boundaries, but shoud not occur also.
cellid=zeros(length(indices),1);
cellid(:)=whichcell;%to keep track to which cell a boundary pixel belongs
totalindices{whichcell,1}= indices ;
totalindicesCELLID{whichcell,1}= cellid ;
overlapimage(indices)=CurrentModelMatrix.ObjectID(whichcell); %we put here on the image the boundaries labeled with cell ID
    end
end
totalindices = cat(1,totalindices{:});
totalindicesCELLID = cat(1,totalindicesCELLID{:});
totalindices(totalindices == 0) = []; 
totalindicesCELLID(totalindicesCELLID == 0) = [];
inside_polygon_ofAllCells = cat(1,PosInsideCell{:});
% 
%  for i = 1:size(IMGCellpatchesFull,1)
%      subplot(4, 4, i)
%      imshow(BoundariesClosed{i})
%  end


indices_with_cellid=[totalindices totalindicesCELLID];
indices_with_cellid=sortrows(indices_with_cellid,1);


analyse_inside_ALLcells=sort(inside_polygon_ofAllCells);


a1=[0;analyse_inside_ALLcells];
b1=[analyse_inside_ALLcells;0];
difflist_area=b1-a1;difflist_area(end) = [];
overlapping_area_indices=analyse_inside_ALLcells(difflist_area==0);


%this is faster
labeled_overlapregions=false(CurrentModelMatrix.rownumber,CurrentModelMatrix.columnnumber);
labeled_overlapregions(overlapping_area_indices)=true;

%toc
    
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Postprocessing. 
%figure;imshow(labeled_overlapregions);impixelinfo
%figure;imshow(label_obarea);impixelinfo
%figure;imshow(area1);impixelinfo
%figure;imshow(area0);impixelinfo
%figure;imshow(image_virtualimage,[]);impixelinfo


%First: remove overlaps which are smaller than 30 pixels
overlapregions_area=regionprops(labeled_overlapregions,'Area');
howmanyoverlaps=length(overlapregions_area);
labeled_overlapregions=bwlabel(labeled_overlapregions);%this is to go faster with regionprops using logical array
for i=1:howmanyoverlaps
    if overlapregions_area(i,1).Area<30
        labeled_overlapregions(labeled_overlapregions==i)=0;
    end
end
labeled_overlapregions=logical(labeled_overlapregions);

%Secondly, remove overlaps in which just boundaries overlap...
overlapregions_area=regionprops(labeled_overlapregions,'Area');
howmanyoverlaps=length(overlapregions_area);
overlapimage_border=overlapimage;
overlapimage_border(overlapimage_border>0)=1;
%figure;imshow(overlapimage_border,[]);impixelinfo
labeled_overlapregions=bwlabel(labeled_overlapregions);
for index_overlap=1:howmanyoverlaps
    ix=overlapimage_border(labeled_overlapregions==index_overlap);
    reducedsize=overlapregions_area(index_overlap,1).Area-length(ix(ix==1)); %ok must check here with the 10000 | (ix>399 & ix<420)
    if reducedsize<2
    labeled_overlapregions(labeled_overlapregions==index_overlap)=0;
    end
end

%this is the greatest trick to replace unique :) It is tremendously faster
%as diff is a built in function
% labeled_overlapregions = bwlabel(labeled_overlapregions);IntraStep = labeled_overlapregions(:);
% howmanyoverlaps = IntraStep([true;diff(labeled_overlapregions(:))>0]);
% howmanyoverlaps = howmanyoverlaps([true;diff(howmanyoverlaps(:))>0]);
labeled_overlapregions=bwlabel(labeled_overlapregions);
howmanyoverlaps=unique(labeled_overlapregions(:));
%
if length(howmanyoverlaps) > 1   
labeled_overlapregions=bwlabel(labeled_overlapregions);
howmanyoverlaps=unique(labeled_overlapregions(:));
howmanyoverlaps(howmanyoverlaps == 0) = [];
howmanyoverlaps = howmanyoverlaps';
 %if it is == to 1, this means that there is no overlap region

%Thirdly, find out with which nuclei they overlap
corresponding_nuclei=zeros(size(howmanyoverlaps,2),2);
count = 0;
pos2Rmv = [];
ObjectstoErase = [];
for i=howmanyoverlaps
    
   ix=find(labeled_overlapregions==i);%indices of pixels forming the overlap region
  
       positionTMP = zeros(size(indices_with_cellid,1),1);
       for iMat = 1:size(ix,1)
       TMP = indices_with_cellid(:,1) == ix(iMat);
       positionTMP = positionTMP + TMP;
       end
       positionTMP = logical(positionTMP);%all positions that correspond to the overlap indices in the indices_with_cellid matrix 
       %now check if this correspons to a single ID
       
       TestOnIDs = indices_with_cellid(positionTMP,2);
%        TestOnIDs = TestOnIDs([true;diff(TestOnIDs(:))>0]);
%        TestOnIDs = TestOnIDs([true;diff(TestOnIDs(:))>0]);             
       TestOnIDs = unique(TestOnIDs); %IDs involved in the overlap.
        
       if length(TestOnIDs)==2  
           if TestOnIDs(1) ~= 0 && TestOnIDs(2) ~= 0
           count = count+1;
       corresponding_nuclei(count,:)=TestOnIDs';
           end
       elseif length(TestOnIDs)==1
           if TestOnIDs(1) ~= 0
       ObjectstoErase = [ObjectstoErase TestOnIDs];
       pos2Rmv = [pos2Rmv i];
           end
       else%basically if there is more than 2
           if TestOnIDs(1) ~= 0 && TestOnIDs(2) ~= 0
           count = count+1;
       corresponding_nuclei(count,:)=TestOnIDs(1:2)'; %treat the case where there is only one oobject detected in the overlap region
           end       
       end
end

for i = ObjectstoErase
 labeled_overlapregions(labeled_overlapregions==i) = 0;   
end
%figure;imshow(labeled_overlapregions,[]);impixelinfo
%figure;imshow(overlapimage,[]);impixelinfo
%not satisfyingly resolved.....

%figure;imshow(labeled_overlapregions,[]);impixelinfo


%Forthly, determine characteristics of the regions...
BWlogical = logical(labeled_overlapregions);
extremalistStruct=regionprops(BWlogical,'Extrema');
centroidOverlapStruct=regionprops(BWlogical,'centroid');    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now determine, which nodes are inside the regions, which nodes overlap
%etc...
%i=9;
%special case, now the question is to know why it inserts 0...
% pos2skip = corresponding_nuclei(:,1) == 0 | corresponding_nuclei(:,2) == 0;
% corresponding_nuclei(pos2skip,:) = []; should not occur anymore
howmanyoverlaps(pos2Rmv) = [];
posToGoTrgh = 1:size(howmanyoverlaps,2);%-sum(pos2skip);%;posToGoTrgh(pos2skip) = [];
%should be fine now but we keep it for the moment
clear howmanyoverlaps
for i=posToGoTrgh
    
    nuclei_involved=corresponding_nuclei(i,:);
    extremalist=extremalistStruct(i,1).Extrema;
    centroidOverlap=centroidOverlapStruct(i,1).Centroid(1,:);
    whichcell=nuclei_involved(1,1);
    checkedcell=nuclei_involved(1,2);
    
    for iCell=[whichcell checkedcell]
    %iCell=whichcell
    %iCell=checkedcell
        if CurrentModelMatrix.Nuclei_Location(iCell,1)-centroidOverlap(1)>0 %the overlapping stuff is shifted above the cell
            if CurrentModelMatrix.Nuclei_Location(iCell,2) - centroidOverlap(1,2)<0 %on the right
               pointlist=[extremalist(1,:); extremalist(8,:);extremalist(4,:);extremalist(5,:)];
            else%on the left or exactly centered on Y
               pointlist=[extremalist(7,:); extremalist(6,:);extremalist(3,:);extremalist(2,:)];
            end
        else%the overlapping stuff is shifted under the cell
            if CurrentModelMatrix.Nuclei_Location(iCell,2) - centroidOverlap(2)<0 %on the right
               pointlist=[extremalist(7,:); extremalist(6,:);extremalist(3,:);extremalist(2,:)];
            else%on the left or exactly centered on Y
               pointlist=[extremalist(1,:); extremalist(8,:);extremalist(4,:);extremalist(5,:)];
            end
        end

        %Calculate the angle in which the "intruding" node is.
        %problem atan only defined on [-pi/2,pi/2]. Hence you must introduce a
        %special treatment of other cases: use artanindex for the cases
        %[-pi/2,pi/2] and arctanindexneg for the other cases [pi/2,3/2pi]
        %finally change the negative angles for the subsequent search (find(max(...)) to positive
        %ones
        angles = nan(size(pointlist,1),1);
        arctanindex= (-CurrentModelMatrix.Nuclei_Location(iCell,1)+pointlist(:,1)) >= 0 ;%find points that are under the centroid 
        arctanindexneg= (-CurrentModelMatrix.Nuclei_Location(iCell,1)+pointlist(:,1)) <0 ;%find points that are over the centroid
        angles(arctanindex)=atan((CurrentModelMatrix.Nuclei_Location(iCell,2)-pointlist(arctanindex,2))./(-CurrentModelMatrix.Nuclei_Location(iCell,1)+pointlist(arctanindex,1)));
        angles(arctanindexneg)=pi- atan((CurrentModelMatrix.Nuclei_Location(iCell,2)-pointlist(arctanindexneg,2))./(CurrentModelMatrix.Nuclei_Location(iCell,1)-pointlist(arctanindexneg,1)));
        angles(angles<0)=2*pi+angles(angles<0); %ok here we got angles
        angles(isnan(angles)) = [];
        %angles=sort(angles); [why did I sort it? => need ordering for assignment of extrema points]

        if CurrentModelMatrix.Nuclei_Location(iCell,1)-centroidOverlap(1)>0 %Over the cell
            if CurrentModelMatrix.Nuclei_Location(iCell,2) - centroidOverlap(2)<0 %on the right
               lowerboundary=min(angles);
               upperboundary=max(angles);
               extremapoint_temp1=pointlist(angles==lowerboundary,:);%could be that unlikely two angles have the same value, let see.
               extremapoint_temp2=pointlist(angles==upperboundary,:);
               extremapoint_temp3=[lowerboundary;upperboundary];
               nodestochange= interactionmatrix.bordernodes{iCell}(:,2)>lowerboundary & interactionmatrix.bordernodes{iCell}(:,2)<upperboundary  ;
               if any(interactionmatrix.bordernodes{iCell}(:,2)<lowerboundary)
                   %find the closest angle that is lower than the lower
                   %boundary of the intruding part of the cell
                   reservenode1= interactionmatrix.bordernodes{iCell}(:,2) == max(interactionmatrix.bordernodes{iCell}(interactionmatrix.bordernodes{iCell}(:,2)<lowerboundary,2));
               else
                   %otherwise we take the last one
                   reservenode1=false(length(interactionmatrix.bordernodes{iCell}(:,2)),1);
                   reservenode1(end) = true;
               end
               
               if any(interactionmatrix.bordernodes{iCell}(:,2)>upperboundary)
                   %find the closest angle that is bigger than the upper
                   %boundary of the intruding part of the cell
                   reservenode2=interactionmatrix.bordernodes{iCell}(:,2)==min(interactionmatrix.bordernodes{iCell}(interactionmatrix.bordernodes{iCell}(:,2)>upperboundary,2));
               else
                   %otherwise we take the first one
                   reservenode2 = false(length(interactionmatrix.bordernodes{iCell}(:,2)),1);
                   reservenode2(1) = true;
               end     
            else %on the left
               lowerboundary=min(angles);
               upperboundary=max(angles);
               extremapoint_temp1=pointlist(angles==lowerboundary,:);
               extremapoint_temp2=pointlist(angles==upperboundary,:);
               extremapoint_temp3=[lowerboundary;upperboundary];
               nodestochange = interactionmatrix.bordernodes{iCell}(:,2)>lowerboundary & interactionmatrix.bordernodes{iCell}(:,2)<upperboundary ;
               
               if any(interactionmatrix.bordernodes{iCell}(:,2)<lowerboundary)
                   reservenode1 = interactionmatrix.bordernodes{iCell}(:,2)==max(interactionmatrix.bordernodes{iCell}(interactionmatrix.bordernodes{iCell}(:,2)<lowerboundary,2));
               else
                   reservenode1 = false(length(interactionmatrix.bordernodes{iCell}(:,2)),1);
                   reservenode1(end) = true;
               end
               
               if any(interactionmatrix.bordernodes{iCell}(:,2)>upperboundary)
                    reservenode2= interactionmatrix.bordernodes{iCell}(:,2)==min(interactionmatrix.bordernodes{iCell}(interactionmatrix.bordernodes{iCell}(:,2)>upperboundary,2));
               else
                   reservenode2 = false(length(interactionmatrix.bordernodes{iCell}(:,2)),1);
                   reservenode2(1) = true;
               end
            end
            else
            if CurrentModelMatrix.Nuclei_Location(iCell,2) - centroidOverlap(2)<0 
              if any(angles<1/2*pi)
                    boundary1=max(angles(angles<1/2*pi));
                    boundary2=min(angles(angles>pi));
                    extremapoint_temp1=pointlist(angles==boundary1,:);
                    extremapoint_temp2=pointlist(angles==boundary2,:);
                    extremapoint_temp3=[boundary1;boundary2];
                    nodestochange = interactionmatrix.bordernodes{iCell}(:,2)>boundary2 | interactionmatrix.bordernodes{iCell}(:,2)<boundary1 ;
                   %nodestochange=find((interactionmatrix.bordernodes{iCell}(:,2)>boundary2)|(interactionmatrix.bordernodes{iCell}(:,2)<boundary1));
                    reservenode1= interactionmatrix.bordernodes{iCell}(:,2)==min(interactionmatrix.bordernodes{iCell}((interactionmatrix.bordernodes{iCell}(:,2)>boundary1),2));
                    reservenode2= interactionmatrix.bordernodes{iCell}(:,2)==max(interactionmatrix.bordernodes{iCell}((interactionmatrix.bordernodes{iCell}(:,2)<boundary2),2));
               else
                    boundary1=min(angles);
                    boundary2=max(angles);
                    extremapoint_temp1=pointlist(angles==boundary1,:);
                    extremapoint_temp2=pointlist(angles==boundary2,:);
                    extremapoint_temp3=[boundary1;boundary2];
                    nodestochange= interactionmatrix.bordernodes{iCell}(:,2)<boundary2 & interactionmatrix.bordernodes{iCell}(:,2)>boundary1 ;
                    if any(interactionmatrix.bordernodes{iCell}(:,2)>boundary2)
                    reservenode2= interactionmatrix.bordernodes{iCell}(:,2)>boundary2;
                    reservenode2= cumsum(reservenode2) == 1 & reservenode2  ; %this gets the firt true of a logical vector, and we do not use find!
                    else
                    reservenode2 = false(length(interactionmatrix.bordernodes{iCell}(:,2)),1);
                    reservenode2(1) = true;
                    end
                    reservenode1= interactionmatrix.bordernodes{iCell}(:,2)>boundary2;
                    reservenode1= cumsum(reservenode1) == max(cumsum(reservenode1)) & reservenode1 ;%this gets the last true of a logical vector, and we do not use find!
                    %reservenode1=find(interactionmatrix.bordernodes{iCell}((interactionmatrix.bordernodes{iCell}(:,2)<boundary1),2), 1, 'last' );
              end
              
            else
              
                if any(angles>3/2*pi)
                    boundary1=min(angles(angles>3/2*pi));
                    boundary2=max(angles(angles<pi));
                    extremapoint_temp1=pointlist(angles==boundary1,:);
                    extremapoint_temp2=pointlist(angles==boundary2,:);
                    extremapoint_temp3=[boundary1;boundary2];
                    nodestochange= interactionmatrix.bordernodes{iCell}(:,2)<boundary2 | interactionmatrix.bordernodes{iCell}(:,2)>boundary1 ;
                    reservenode1=interactionmatrix.bordernodes{iCell}(:,2) == min(interactionmatrix.bordernodes{iCell}((interactionmatrix.bordernodes{iCell}(:,2)>boundary2),2));
                    reservenode2=interactionmatrix.bordernodes{iCell}(:,2) == max(interactionmatrix.bordernodes{iCell}((interactionmatrix.bordernodes{iCell}(:,2)<boundary1),2));
                else
                    boundary1=min(angles);
                    boundary2=max(angles(angles<pi));
                    extremapoint_temp1=pointlist(angles==boundary1,:);
                    extremapoint_temp2=pointlist(angles==boundary2,:);
                    extremapoint_temp3=[boundary1;boundary2];
                    nodestochange= interactionmatrix.bordernodes{iCell}(:,2) < boundary2 & interactionmatrix.bordernodes{iCell}(:,2) > boundary1 ;
                    if any(interactionmatrix.bordernodes{iCell}(:,2)<boundary1)
                    reservenode1 = interactionmatrix.bordernodes{iCell}(:,2)<boundary1 ;
                    reservenode1= cumsum(reservenode1) == max(cumsum(reservenode1)) & reservenode1;                
                    else
                        reservenode1=false(length(interactionmatrix.bordernodes{iCell}(:,2)),1);
                        reservenode1(end) = true;
                    end
                    reservenode2 = interactionmatrix.bordernodes{iCell}(:,2)>boundary2;
                    reservenode2= cumsum(reservenode2) == 1 & reservenode2 ;
                end 
            end
        end
        interactionmatrix.nodestochange{iCell}(nodestochange)=1;
        if iCell==whichcell
        interactionmatrix.bordernodes{iCell,1}(nodestochange,7)=checkedcell;
        else
        interactionmatrix.bordernodes{iCell,1}(nodestochange,7)=whichcell;
        end
        
        %let's restrict to the closest reserve node
        if sum(reservenode1)>1 && any(nodestochange)
            %we look at the first position
            NTCcomp = cumsum(nodestochange) == 1 & nodestochange;
            NTCref = 1:size(NTCcomp,1); RNlist = NTCref(reservenode1);NTCpos = NTCref(NTCcomp);
            %
            reservenode1 = abs(RNlist(:)-NTCpos) == min(abs(RNlist(:)-NTCpos));
            NTCref(:) = false;
            NTCref(RNlist(reservenode1)) = true; reservenode1 = NTCref';
        end
        if sum(reservenode2) > 1 && any(nodestochange)
            NTCcomp = cumsum(nodestochange) == 1  & nodestochange ;%keep the first one
            NTCref = 1:size(NTCcomp,1); RNlist = NTCref(reservenode2);NTCpos = NTCref(NTCcomp);%find  positions
            %
            reservenode2 = abs(RNlist(:)-NTCpos) == min(abs(RNlist(:)-NTCpos));%defines the closest positions
            NTCref(RNlist(reservenode2)) = true; reservenode2 = NTCref';
            %reservenode2=reservenode2(abs(reservenode2(:)-nodestochange(1))==min(abs(reservenode2(:)-nodestochange(1))));%keep
            %as areminder of what was here before
        end
        
        if sum(reservenode1)>1
            reservenode1 = cumsum(reservenode1) == 2 & reservenode1;%double cumsum allows to stay away from find and avoid problems of streches of 0s in logical arrays
        end
        if sum(reservenode2)>1 
            reservenode2 = cumsum(reservenode2) == 2 & reservenode2;
        end
        
%         if isempty(reservenode1) 
%             reservenode1=[];%in the mean time to have better fix
%         end
%         if isempty(reservenode2)
%             reservenode2=[];
%         end
        
        if iCell==whichcell
            interactionmatrix.reservenode{iCell,checkedcell}=[reservenode1,reservenode2];
            if ~isempty(interactionmatrix.extremaList{iCell,checkedcell})
            templength=length(interactionmatrix.extremaList{iCell,checkedcell}(:,1));
            else
            templength=0;
            end
            interactionmatrix.extremaList{iCell,checkedcell}(1+templength,1:2)=extremapoint_temp1(1,:);
            interactionmatrix.extremaList{iCell,checkedcell}(2+templength,1:2)=extremapoint_temp2(1,:);
            interactionmatrix.extremaList{iCell,checkedcell}(templength+1:templength+2,3)=extremapoint_temp3;
        else
            interactionmatrix.reservenode{iCell,whichcell}=[reservenode1,reservenode2];
            if ~isempty(interactionmatrix.extremaList{iCell,whichcell})
            templength=length(interactionmatrix.extremaList{iCell,whichcell}(:,1));
            else
            templength=0;
            end
        interactionmatrix.extremaList{iCell,whichcell}(1+templength,1:2)=extremapoint_temp1(1,:);
        interactionmatrix.extremaList{iCell,whichcell}(2+templength,1:2)=extremapoint_temp2(1,:);
        interactionmatrix.extremaList{iCell,whichcell}(templength+1:templength+2,3)=extremapoint_temp3;
        end
        
    end %over pair of cells
      
end %over all areas
end


temp4=cell(CurrentModelMatrix.numberofcells,1);
for i=1:CurrentModelMatrix.numberofcells
    
    test = interactionmatrix.nodestochange{i}==1;
    if any(test)
    temp=sum(test);
    temp2=zeros(temp,1);
    temp2(1:end,1)=i;
    temp3=[find(test), temp2];
    temp4{i}=temp3;
    end
end
temp4 = cat(1,temp4{:});

%virtualimage=visualize_celloutline_model_cell(interactionmatrix,CurrentModelMatrix);
%plot_overlapping_nodes_reserve(CurrentModelMatrix, virtualimage, interactionmatrix, temp4)




%######comment 1

%####end






% #Comment 2
% %so now look in the bounding box of the cell
% not_inside_total1=nan(20.*((box_coord(iCell,2) - box_coord(iCell,1))+1),1); %we assume 20 is enough allocation
% posIX = 0;
% for i=box_coord(iCell,1):box_coord(iCell,2)
%   
%     ix_coor=find(overlapimage(box_coord(iCell,3):box_coord(iCell,4),i)==CurrentModelMatrix.ObjectID(iCell,1));%do we find the actual cell boundaries in our piece of image
%     ix_8000=find(overlapimage(box_coord(iCell,3):box_coord(iCell,4),i)==8000);% do we find overlapping pixels
%     ix_8000coord=ix_8000+(box_coord(iCell,3)-1)+(i-1)*CurrentModelMatrix.rownumber;%absolute position like given using A=FIND(B).    
%     ix_coor2 = nan(length(ix_8000),1);%how many overlaping pixels?
%     count = 0;
%     for index_n=1:length(ix_8000)%for all overlapping pixels
%     ix_8000coord_transform=indices_with_cellid((ix_8000coord(index_n)==indices_with_cellid(:,1)),2);%finds back the IDs of cells that are overlapping
%     if find(iCell==ix_8000coord_transform)%if you find pixels belonging to the cell you are looking at
%         count = count+1;
%         ix_coor2(count)= ix_8000(index_n);%if the overlap concerns our cell store it.
%     end
%     end
%     ix_coor2(isnan(ix_coor2)) = [];
%     
%     
%     ix=[ix_coor;ix_coor2];%ix_coor : boundaries of our current cell, ix_coor2: boundaries of our current cell that overlap with other cells
%     ix=sort(ix);%logically, the ovelapping pixels get cited 2 times 
%     %general coordinates
%     ix1=ix(1)+(box_coord(iCell,3)-1);
%     ix2=ix(length(ix))+(box_coord(iCell,3)-1);
%     %
%     not_inside1=box_coord(iCell,3):(ix1-1);
%     not_inside2=(ix2+1):box_coord(iCell,4);
%     
%     %absolute coordinates
%     not_inside=[not_inside1 not_inside2];
%     not_inside=CurrentModelMatrix.rownumber*(i-1)+not_inside;
%     %
%     
%     bdrTMp = size(not_inside,2);
%     not_inside_total1(posIX+1:posIX+bdrTMp,1)= not_inside';
%     posIX = posIX + bdrTMp;
%    
% end
% not_inside_total1(isnan(not_inside_total1)) = [];
% 
% 
% not_inside_total2= nan(20.*((box_coord(iCell,4) - box_coord(iCell,3))+1),1);
% posIX = 0;
% %same logic here
% for i=box_coord(iCell,3):box_coord(iCell,4)
%     
%     ix_coor=find(overlapimage(i,(box_coord(iCell,1):box_coord(iCell,2)))==CurrentModelMatrix.ObjectID(iCell,1));
%     ix_8000=find(overlapimage(i,box_coord(iCell,1):box_coord(iCell,2))==8000);
%     ix_8000coord=(((ix_8000-1)+box_coord(iCell,1))-1)*CurrentModelMatrix.rownumber+i;
%     ix_coor2 = nan(length(ix_8000),1);
%     count = 0;
%     for index_n=1:length(ix_8000)
%     ix_8000coord_transform=indices_with_cellid((ix_8000coord(index_n)==indices_with_cellid(:,1)),2);
%     if find(iCell==ix_8000coord_transform)
%         count = count+1;
%         ix_coor2(count)= ix_8000(index_n);
%     end
%     end
%     ix_coor2(isnan(ix_coor2)) = [];
%     
%     ix=[ix_coor';ix_coor2];
%     ix=sort(ix);
%     ix1=ix(1)+box_coord(iCell,1)-1;
%     ix2=ix(length(ix))+box_coord(iCell,1)-1;
%     not_inside1=box_coord(iCell,1):(ix1-1);
%     not_inside2=(ix2+1):box_coord(iCell,2);
%     not_inside=[not_inside1,not_inside2];
%     not_inside=CurrentModelMatrix.rownumber*(not_inside-1)+i;
%     bdrTMp = size(not_inside,2);
%     not_inside_total2(posIX+1:posIX+bdrTMp,1)= not_inside';
%     posIX = posIX+bdrTMp;
%     
% end
% not_inside_total2(isnan(not_inside_total2)) = [];
% 
% %what follows can be done in a faster way
% all_pixels = nan(20.*((box_coord(iCell,2) - box_coord(iCell,1))+1),1);
% posIX=0;
% for i=box_coord(iCell,1):box_coord(iCell,2)
%     
%     index_list=box_coord(iCell,3):box_coord(iCell,4);
%     index_list=(i-1)*CurrentModelMatrix.rownumber+index_list;
%     bdrTMp = size(index_list,2);
%     all_pixels(posIX+1:posIX+bdrTMp,1)= index_list';
%     posIX = posIX + bdrTMp;
%     
% end
% all_pixels(isnan(all_pixels)) = [];
% 
% not_inside_total=unique([not_inside_total1;not_inside_total2]);%this list all positions that are not in our object iCell
% % overlaptest=overlapimage;
% inside_polygon=unique(setdiff(all_pixels,not_inside_total));%list pixels inside the object
% 
% inside_polygon_ofAllCells=[inside_polygon_ofAllCells;inside_polygon];%for all cells
% end
% end


%overlaptest=overlapimage;

%overlaptest(inside_polygon_total)=400;
%figure;imshow(overlaptest,[]);impixelinfo
%EndComment#2  


%  Node2LookAtXRel = Node2LookAtX - (box_coord(whichcell,1)-1);%minus 1 of course as a coordinate equal to 0 is a nonsens
%  Node2LookAtYRel = Node2LookAtY - (box_coord(whichcell,3)-1);
%  Node2LookAtXRel(Node2LookAtXRel <= 0) = 1;
%  Node2LookAtYRel(Node2LookAtYRel <= 0) = 1;
%  TMPimg = false(max(Node2LookAtXRel),max(Node2LookAtYRel));
%  TMPpos = ((Node2LookAtYRel-1) .* max(Node2LookAtXRel)) + Node2LookAtXRel;
%  TMPimg(TMPpos) = true;
%  TMPimg = bwconvhull(TMPimg);
%  TMPimg = bwperim(TMPimg);
%  [X Y] = find(TMPimg);
%  clear TMPimg TMPpos
%  %replace the object with real coordinates
%  X = X+(box_coord(whichcell,1)-1);
%  Y = Y+(box_coord(whichcell,3)-1);
%  X(X > CurrentModelMatrix.rownumber) = CurrentModelMatrix.rownumber;
%  Y(Y > CurrentModelMatrix.columnnumber) = CurrentModelMatrix.columnnumber;
%  indices = ((Y-1).*CurrentModelMatrix.rownumber) + X;
%  %to see what was before go to the end comment 1.


% for i=1:size(ParForList,1)
%    inde = ParForList(i,1);   
%     if inde==lastVal %for the last node in the list
%           newindices=linept_FAST(CurrentModelMatrix,Node2LookAtX(inde), ... 
%           Node2LookAtY(inde),Node2LookAtX(FirstVal),Node2LookAtY(FirstVal));
%           indices{i} = newindices ;
% %         bndrTMP = size(newindices,1);
% %         indices(currentIX+1:currentIX+bndrTMP,1) = newindices ;
% %         currentIX = bndrTMP + currentIX;            
%     else %for the others
%           indePlus1 = ParForList(i+1,1);
%           newindices=linept_FAST(CurrentModelMatrix,Node2LookAtX(inde), ... 
%           Node2LookAtY(inde),Node2LookAtX(indePlus1),Node2LookAtY(indePlus1));
%           indices{i} = newindices ;
% %         bndrTMP = size(newindices,1);        
% %         indices(currentIX+1:currentIX+bndrTMP,1) = newindices ; %indices of points to reconstitute boundaries of cells
% %         currentIX = bndrTMP + currentIX;
%     end
% end

% indices = cat(1,indices{:});


% a=[0;indices_with_cellid(:,1)];
% b=[indices_with_cellid(:,1);0];
% difflist=b-a;
% overlapping_line_indices= difflist==0 ;
% overlapping_line_points=indices_with_cellid(overlapping_line_indices,1);
% 
% overlapimage(overlapping_line_points)=8000; %put more intense signal where it overlaps
% %figure;imshow(overlapimage,[]);impixelinfo;
% %%%%%




% %calculation
% for iCell=1:CurrentModelMatrix.numberofcells
%  if ~CurrentModelMatrix.FreezeTag(iCell)
%     ID = CurrentModelMatrix.ObjectID(iCell);
% IMGCellpatches{iCell} = overlapimage(box_coord(iCell,1):box_coord(iCell,2),box_coord(iCell,3):box_coord(iCell,4));
% %this matrix is store at the very beggining during init. step
% IdOverlapPix = AbsolutePos{iCell}(IMGCellpatches{iCell} == 8000); 
% PosTMP = arrayfun(@(x) x == indices_with_cellid(:,1),IdOverlapPix,'UniformOutput',false);PosTMP = cat(2,PosTMP{:});PosTMP = max(PosTMP,[],2);
% ObjectConcernedByOL = indices_with_cellid(PosTMP,:); %here we get cells implicated in OL.
% IdOverlapPix4iCell{iCell} = ObjectConcernedByOL(ObjectConcernedByOL(:,2) == iCell,1);%I filterOut just the ones implicated in OL with iCell
% BoundariestoFill{iCell} = IMGCellpatches{iCell} == ID ;
% for iClose = 1:size(IdOverlapPix4iCell{iCell},1)
%     BoundariestoFill{iCell}(AbsolutePos{iCell} == IdOverlapPix4iCell{iCell}(iClose)) = true;
% end
% %important to close the object even where there is an overlap involving iCell
% BoundariesClosed{iCell} = bwmorph(BoundariestoFill{iCell},'bridge');
% IMGCellpatchesFull{iCell} = imfill(BoundariesClosed{iCell},'holes');
% % = bwconvhull(IMGtmp);
% tostore = AbsolutePos{iCell}(IMGCellpatchesFull{iCell});
% PosInsideCell{iCell,1} = tostore;
%  end
% end
% inside_polygon_ofAllCells = cat(1,PosInsideCell{:}); 
% %%%%to see what was before just go at the end comment#2
% for i = 1:size(IMGCellpatchesFull,1)
%     subplot(4, 4, i)
%     imshow(IMGCellpatchesFull{i})
% end
