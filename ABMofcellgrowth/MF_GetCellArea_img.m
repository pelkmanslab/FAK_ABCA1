function[totalarea]=MF_GetCellArea_img(CurrentModelMatrix,interactionmatrix,labels)

%[MF]Totally changed 09/07/13

%gives back an image with the area occupied by cells from the set of nodes

%let's store this for later
PosInsideCell = cell(CurrentModelMatrix.numberofcells,1);%positions of the current image that are in the object concerned
AbsolutePos = cell(CurrentModelMatrix.numberofcells,1);%the actual patch of the object of interrest, with absolute pos in the total image > fast indexing
IMGCellpatchesFull = cell(CurrentModelMatrix.numberofcells,1);
BoundariesClosed = cell(CurrentModelMatrix.numberofcells,1);
% IdOverlapPix4iCell = cell(CurrentModelMatrix.numberofcells,1);%same patch always, but with only overlapping pixels involving the current object of interrest
% IMGCellpatches = cell(CurrentModelMatrix.numberofcells,1);%the actual patch of the object of interrest, contains the object boundaries as objID pixels 
% BoundariestoFill = cell(CurrentModelMatrix.numberofcells,1);%same patch always, with outlines in a logical array
totalarea=zeros(CurrentModelMatrix.rownumber,CurrentModelMatrix.columnnumber);%hum...
nodesofcell=interactionmatrix.bordernodes;

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


leftboundary=min(interactionmatrix.bordernodes{whichcell,1}(node_list{whichcell,1},3));
rightboundary=max(interactionmatrix.bordernodes{whichcell,1}(node_list{whichcell,1},3));
topboundary=min(interactionmatrix.bordernodes{whichcell,1}(node_list{whichcell,1},4));
lowerboundary=max(interactionmatrix.bordernodes{whichcell,1}(node_list{whichcell,1},4));
box_coord(whichcell,:)=[topboundary,lowerboundary,leftboundary,rightboundary];

%let draw lines between nodes...
indices = cell(size(node_list{whichcell,1},1),1);
Node2LookAtX = nodesofcell{whichcell}(node_list{whichcell,1},4);
Node2LookAtY = nodesofcell{whichcell}(node_list{whichcell,1},3);
ParForList = (1:size(Node2LookAtX,1))' ;
FirstVal = ParForList(1,1);
lastVal = ParForList(end,1);

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
PosInsideCell{whichcell} = AbsolutePos{whichcell}(IMGCellpatchesFull{whichcell});

    end
end

if labels == false
inside_polygon_ofAllCells = cat(1,PosInsideCell{:});
totalarea(inside_polygon_ofAllCells)=true;
elseif labels == true
 for i = 1:size(PosInsideCell,1)
 totalarea(PosInsideCell{i}) =  CurrentModelMatrix.ObjectID(i);  
 end       
end
  