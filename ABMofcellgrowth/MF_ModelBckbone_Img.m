function [virtualimage]=MF_ModelBckbone_Img(interactionmatrix,CurrentModelMatrix)

%%%%%%%%%%%%%%
%visualize points around cell nucleus

%draw lines between points but only those lines where the nodes are in the
%image and not outside the boundaries (nodes can be placed outside by
%randomly assigning nodes around a nucleus centre. Hence those cases must
%be excluded from drawing since otherwise the image dimensions of
%virtualimage will change.
virtualimage=MF_GetNucleiFromVectData_img(CurrentModelMatrix);
% PosInsideCell = cell(CurrentModelMatrix.numberofcells,1);%positions of the current image that are in the object concerned
AbsolutePos = cell(CurrentModelMatrix.numberofcells,1);%the actual patch of the object of interrest, with absolute pos in the total image > fast indexing
IMGCellpatchesFull = cell(CurrentModelMatrix.numberofcells,1);
BoundariesClosed = cell(CurrentModelMatrix.numberofcells,1);
nodesofcell=interactionmatrix.bordernodes;

totalindices = cell(CurrentModelMatrix.numberofcells,1);
box_coord=zeros(CurrentModelMatrix.numberofcells,4);
node_list = cell(CurrentModelMatrix.numberofcells,1);
CoorXobj = cell(CurrentModelMatrix.numberofcells,1);
CoorYobj = cell(CurrentModelMatrix.numberofcells,1);
CoorYRel = cell(CurrentModelMatrix.numberofcells,1);
CoorXRel = cell(CurrentModelMatrix.numberofcells,1);



for whichcell=1:CurrentModelMatrix.numberofcells
    node_list{whichcell,1}=false(interactionmatrix.numberofbordernodes(whichcell),1);
%[MF] check for the eligibility of the node, checks if in image in
%particular.
for inde=1:interactionmatrix.numberofbordernodes(whichcell)
   if (nodesofcell{whichcell}(inde,3)>=1 && nodesofcell{whichcell}(inde,4)>=1 ...
       && nodesofcell{whichcell}(inde,4)<= CurrentModelMatrix.rownumber ... 
       && nodesofcell{whichcell}(inde,3)<= CurrentModelMatrix.columnnumber)
   %if the node is in the image we can process it, if not freeze the cell.
   %we can come on this later
   node_list{whichcell,1}(inde,1)=true; %if the node can be treated, it gets a number
   end
end


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
BoundariesClosed{whichcell} = bwmorph(BoundariesClosed{whichcell},'spur');
% BoundariesClosed{whichcell} = bwmorph(BoundariesClosed{whichcell},'remove');
% IMGCellpatchesFull{whichcell} = imfill(BoundariesClosed{whichcell},'holes');

indices = AbsolutePos{whichcell}(BoundariesClosed{whichcell});
% PosInsideCell{whichcell} = AbsolutePos{whichcell}(IMGCellpatchesFull{whichcell});
indices(indices == 0) = []; %actually this should not occur.
totalindices{whichcell,1}= indices ;

end
pos = cat(1,totalindices{:});

virtualimage(pos)=8000;

%plot nodes. Also here nodes which are outside the image dimensions have to
%be excluded for plotting reasons.

%[r,c]=find(interactionmatrix.nodestochange==1);
%c=CurrentModelMatrix.ObjectID(c);
%[r,c]
% for iCell = 1: CurrentModelMatrix.numberofcells
% for i = 1:interactionmatrix.numberofbordernodes(iCell)
% coordinate_condition3=nodesofcell{iCell}(i,3);
% coordinate_condition4=nodesofcell{iCell}(i,4);
% checkcoordinate_condition=~reshape((coordinate_condition3 > 1 & coordinate_condition3 < CurrentModelMatrix.columnnumber & coordinate_condition4 > 1 & coordinate_condition4 < CurrentModelMatrix.rownumber),[],1);
% coordinatesfor1=reshape(CurrentModelMatrix.rownumber.*nodesofcell{iCell}(i,3)+nodesofcell{iCell}(i,4),[],1);
% coordinatesfor1(checkcoordinate_condition)=[];
% virtualimage(coordinatesfor1)=8000+i;
% end
% end
% clear coordinatesfor1 checkcoordinate_condition coordinate_condition3 coordinate_condition4

%testimage =zeros(CurrentModelMatrix.rownumber, CurrentModelMatrix.columnnumber);
%for in=1:CurrentModelMatrix.numberofcells
%testimage= testimage+poly2mask(nodesofcell(:,3,in),nodesofcell(:,4,in),CurrentModelMatrix.rownumber,CurrentModelMatrix.columnnumber);
%end
%figure
%imshow(testimage,[])
%figure; imshow(virtualimage,[]);impixelinfo

