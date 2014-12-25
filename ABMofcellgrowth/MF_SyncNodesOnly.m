function[interactionmatrix,CurrentModelMatrix,energymatrix, output4]=MF_SyncNodesOnly(MovementPar,energymatrix,CurrentModelMatrix,interactionmatrix, overlappinglist,nodeGrowthlimit)
%[MF]shortly checked 06/13
%overlappinglist=overlaplist;
LengthOfNode = floor(nodeGrowthlimit * CurrentModelMatrix.growthtimestepincrement);
minimal_dist_to_nucleimembrane=3;
output4=0;
[tempindex,~]=size(overlappinglist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if there is an overlap, first check with the corresponding node, whether
%it is an innernode, or if the node is already at minimal distance...

%safety step
for i=1:tempindex
iCell=overlappinglist(i,2);
nbNode=overlappinglist(i,1);
%check and erase
if size(interactionmatrix.bordernodes{iCell},1) < nbNode
   overlappinglist(i,:) = []; 
end
end

indextodelete=[];
for i=1:tempindex
iCell=overlappinglist(i,2);
nbNode=overlappinglist(i,1);
if interactionmatrix.bordernodes{iCell,1}(nbNode,5)>0
  indextodelete=[indextodelete,i];
end
end
overlappinglist(indextodelete,:)=[];

[tempindex,~]=size(overlappinglist);

%movecells if required
MeanDistrib = MovementPar(1,1);
StdDistrib = MovementPar(1,2);
movelefttoright=round(random('normal',MeanDistrib,StdDistrib));
movetopdown=round(random('normal',MeanDistrib,StdDistrib));
for i=1:tempindex
iCell=overlappinglist(i,2);
CurrentModelMatrix.Nuclei_Location(iCell,:)=CurrentModelMatrix.Nuclei_Location(iCell,:)+[movelefttoright,movetopdown];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now move all nodes which can be moved

if ~isempty(overlappinglist)
    output4=1;
%  CurrentModelMatrix.nuclei_area=Display_nuclei_ellipsewithrotation_nodisplay(CurrentModelMatrix);
[nuclei_outlines,CurrentModelMatrix.nuclei_area]=MF_GetNucleiOutlinesFromVectData_img(CurrentModelMatrix);
occupied_area=MF_GetCellArea_img(CurrentModelMatrix,interactionmatrix,true);
%  imshow(nuclei_outlines)
%figure;imshow(occupied_area,[]);impixelinfo

%update node position, so that all nodes can be deleted which are not
%"moving nodes"
nbNodeList=zeros(tempindex,1);
for i=1:tempindex
nbNode_temp=overlappinglist(i,1);   
nbNodeList(i,1)=length(find(interactionmatrix.bordernodes{overlappinglist(i,2)}(1:nbNode_temp,5)==0));
end

%Place here all nodes newly and update the interactionmatrix (indCase=3)
for indCase=1:tempindex
     
iCell=overlappinglist(indCase,2);
nbNode=nbNodeList(indCase,1);

interactionmatrix.bordernodes{iCell}(interactionmatrix.bordernodes{iCell}(:,5)>0,:)=[];
interactionmatrix.numberofbordernodes(iCell,1)=interactionmatrix.numberofnodes;

    if ~CurrentModelMatrix.FreezeTag(iCell)
%potential: model that a cell keeps its original direction:calculated per
%cell
cell_vector = [movelefttoright,movetopdown];
energymatrix.current_localPt=MF_SetNodePotential(1,LengthOfNode);
energymatrix.current_movementPt=MF_SetMovementPotential(cell_vector,1,LengthOfNode);
%potential as distance from nucleus and other cells: calculated per cell
[energymatrix]=MF_SetWholeCellPotential(CurrentModelMatrix, interactionmatrix, iCell,occupied_area, CurrentModelMatrix.nuclei_area,nuclei_outlines,energymatrix,1,false);
[resulting_potential, interactionmatrix]=MF_MergePotentials(CurrentModelMatrix, interactionmatrix,energymatrix, iCell,nbNode,LengthOfNode);
%choose random place below the value of the potential function at this
%place...
if ~exist('energymatrix.HasitJustDivided','var')
energymatrix.HasitJustDivided = false(CurrentModelMatrix.numberofcells,1);
end

if energymatrix.HasitJustDivided(iCell) == false
potential_placesList=find(resulting_potential<energymatrix.nodeenergy{iCell}(nbNode,1));
else
potential_placesList=find(resulting_potential<0.9);    
end

if ~isempty(potential_placesList)
    size_area_choosefrom=length(potential_placesList);
    chosen_new_place_prob=round(random('Uniform',0.5,size_area_choosefrom+0.5));
    chosen_new_place=potential_placesList(chosen_new_place_prob);
    [r,~]=size(resulting_potential);
    coord1=mod(chosen_new_place,r);
    if coord1==0
    coord1=1;%rownumber
    coord2=floor(chosen_new_place/r);%columnnumber
    else
    coord2=floor(chosen_new_place/r)+1;
    end
else
    coord1=16;
    coord2=16;
end

%update the interactionmatrix etc...according to the new position...
if LengthOfNode < 15
relative_to_old1=coord1-(LengthOfNode+1);%up-down direction
relative_to_old2=coord2-(LengthOfNode+1);%left-right direction
else
relative_to_old1=coord1-16;%up-down direction
relative_to_old2=coord2-16;%left-right direction
end

interactionmatrix.bordernodes{iCell,1}(nbNode,3)=interactionmatrix.bordernodes{iCell,1}(nbNode,3)+relative_to_old2;
interactionmatrix.bordernodes{iCell,1}(nbNode,4)=interactionmatrix.bordernodes{iCell,1}(nbNode,4)+relative_to_old1;


coordinate=[];
coordinate(:,2)=interactionmatrix.bordernodes{iCell,1}(:,3);
coordinate(:,1)=interactionmatrix.bordernodes{iCell,1}(:,4);
%update now interactionmatrix for all nodes:
interactionmatrix.bordernodes{iCell,1}(:,1)=sqrt((CurrentModelMatrix.Nuclei_Location(iCell,1)-coordinate(:,2)).^2+(CurrentModelMatrix.Nuclei_Location(iCell,2)-coordinate(:,1)).^2);
%calculate angle
angles=zeros(interactionmatrix.numberofnodes,1);
arctanindex=find(-CurrentModelMatrix.Nuclei_Location(iCell,1)+coordinate(:,2)>0);
arctanindexneg=find(-CurrentModelMatrix.Nuclei_Location(iCell,1)+coordinate(:,2)<0);
arctanindexzero=find(-CurrentModelMatrix.Nuclei_Location(iCell,1)+coordinate(:,2)==0);
angles(arctanindex)=atan((CurrentModelMatrix.Nuclei_Location(iCell,2)-coordinate(arctanindex,1))./(-CurrentModelMatrix.Nuclei_Location(iCell,1)+coordinate(arctanindex,2)));
angles(arctanindexneg)=pi- atan((CurrentModelMatrix.Nuclei_Location(iCell,2)-coordinate(arctanindexneg,1))./(CurrentModelMatrix.Nuclei_Location(iCell,1)-coordinate(arctanindexneg,2)));

for index_arctanindex=1:length(arctanindexzero)
    if -CurrentModelMatrix.Nuclei_Location(iCell,2)+coordinate(arctanindexzero(index_arctanindex),1)>0
    angles(arctanindexzero(index_arctanindex))= pi * 3/2;
    else
    angles(arctanindexzero(index_arctanindex)) = pi/2;
    end
end
angles(angles<0)=2*pi+angles(angles<0);
interactionmatrix.bordernodes{iCell}(:,2)=angles;
%calculate here the minimal distance 
factorl=zeros(interactionmatrix.numberofnodes,1);
dist1=factorl;
dist2=factorl;

for node=1:interactionmatrix.numberofnodes
    factort=atan(CurrentModelMatrix.radius(iCell,1)/CurrentModelMatrix.radius(iCell,2)*tan(angles(node)-CurrentModelMatrix.angle(iCell,1)));
    factorl(node,1)=sqrt((CurrentModelMatrix.radius(iCell,1).*cos(factort)).^2 + (CurrentModelMatrix.radius(iCell,2).*sin(factort)).^2);
    dist1(node,1)=factorl(node)*cos(angles(node));
    dist2(node,1)=factorl(node)*sin(angles(node));
    interactionmatrix.bordernodes{iCell}(node,6)=factorl(node,1)+minimal_dist_to_nucleimembrane;
    %interactionmatrix.centertonucleimembrane(node,1,i)=factorl(node,1);
end
    else
[interactionmatrix] = MF_ResetNodePosition(interactionmatrix,CurrentModelMatrix, iCell,1);
    end
end %over list

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now you have to place inner nodes and sort the arrays.
CellListReplaced=unique(overlappinglist(:,2))';
for iCell=CellListReplaced   
%sort rows
interactionmatrix.bordernodes{iCell}=sortrows(interactionmatrix.bordernodes{iCell},2);
[~,interactionmatrix] = MF_BendMembrane(interactionmatrix,CurrentModelMatrix, iCell);
interactionmatrix.bordernodes{iCell}=sortrows(interactionmatrix.bordernodes{iCell},2);  
end

disp('some nodes newly placed in move_overlapping_nodes_again_FAST_STD')
end






