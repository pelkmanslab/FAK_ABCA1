function[occupied_area]=MF_GetCellAreaAroundSelectedNuclei_img(CurrentModelMatrix,interactionmatrix, whichcell,labels)

%[MF]massive changes 06/13: tough is not enough.
%gives back an image with the area occupied by the cells as ones


%%%%%%%%%check which cells are in the vicinity of the given cell%%%%%%%%%%%
%span a window of 300x300 pixels and find all centres of
%nuclei within them

leftright=find((CurrentModelMatrix.Nuclei_Location(:,1)>(CurrentModelMatrix.Nuclei_Location(whichcell,1)-150)) & (CurrentModelMatrix.Nuclei_Location(:,1)<CurrentModelMatrix.Nuclei_Location(whichcell,1)+150));
updown=find((CurrentModelMatrix.Nuclei_Location(:,2)>(CurrentModelMatrix.Nuclei_Location(whichcell,2)-150)) & (CurrentModelMatrix.Nuclei_Location(:,2)<CurrentModelMatrix.Nuclei_Location(whichcell,2)+150));
nucleiinrange=intersect(leftright,updown);
nucleiinrange=setdiff(nucleiinrange, whichcell);
nucleiinrange=nucleiinrange';
occupied_area=zeros(CurrentModelMatrix.rownumber,CurrentModelMatrix.columnnumber);%hum...
clear leftright updown

%let's store this for later
PosInsideCell = cell(length(nucleiinrange),1);%positions of the current image that are in the object concerned
AbsolutePos = cell(length(nucleiinrange),1);%the actual patch of the object of interrest, with absolute pos in the total image > fast indexing
IMGCellpatchesFull = cell(length(nucleiinrange),1);
BoundariesClosed = cell(length(nucleiinrange),1);
CoorXobj = cell(CurrentModelMatrix.numberofcells,1);
CoorYobj = cell(CurrentModelMatrix.numberofcells,1);
CoorYRel = cell(CurrentModelMatrix.numberofcells,1);
CoorXRel = cell(CurrentModelMatrix.numberofcells,1);


if ~isempty(nucleiinrange)
%make an image with all boundaries
nodesofcell=interactionmatrix.bordernodes;
%make an image with all boundaries
% overlapimage=occupied_area;
%totalindices = zeros(560.*CurrentModelMatrix.numberofcells,1);%7 nodes * 80 as maxi line length
% totalindices = cell(CurrentModelMatrix.numberofcells,1);
% totalindicesCELLID = cell(CurrentModelMatrix.numberofcells,1);
box_coord=zeros(CurrentModelMatrix.numberofcells,4);
node_list = cell(CurrentModelMatrix.numberofcells,1);


for whichcell=nucleiinrange
    node_list{whichcell,1}=false(interactionmatrix.numberofbordernodes(whichcell),1);
    if ~CurrentModelMatrix.FreezeTag(whichcell)
%[MF] check for the eligibility of the node, checks if in image in
%particular.

for inde=1:interactionmatrix.numberofbordernodes(whichcell)
   if (nodesofcell{whichcell}(inde,3)>=1 && nodesofcell{whichcell}(inde,4)>=1 ...
       && nodesofcell{whichcell}(inde,4)<= CurrentModelMatrix.rownumber... 
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

% indices = AbsolutePos{whichcell}(BoundariesClosed{whichcell});
PosInsideCell{whichcell} = AbsolutePos{whichcell}(IMGCellpatchesFull{whichcell});
% indices(indices == 0) = []; %actually this should not occur.
    end
end
% totalindices = cat(1,totalindices{:});
% totalindicesCELLID = cat(1,totalindicesCELLID{:});
% totalindices(totalindices == 0) = []; 
% totalindicesCELLID(totalindicesCELLID == 0) = [];
if labels == false
inside_polygon_ofAllCells = cat(1,PosInsideCell{:});
occupied_area(inside_polygon_ofAllCells)=true;
elseif labels == true
 for i = 1:size(PosInsideCell,1)
 occupied_area(PosInsideCell{i}) =  CurrentModelMatrix.ObjectID(i);  
 end       
end
end   

