function[interactionmatrix,CurrentModelMatrix,energymatrix]=MF_SyncNodesWithRandNuclWiggling(energymatrix,CurrentModelMatrix,interactionmatrix,nodeGrowthlimit)
%[MF] I must merge it with model_move_nodes_STD
minimal_dist_to_nucleimembrane=3;
LengthOfNode = floor(nodeGrowthlimit * CurrentModelMatrix.growthtimestepincrement);
ix=[];
for index1=1:CurrentModelMatrix.numberofcells
    for index2=1:CurrentModelMatrix.numberofcells
        ix1=[];
if~isempty(interactionmatrix.extremaList{index1,index2})
ix1=index1;
end
ix=[ix,ix1];
    end
end

%potential: model conservativeness around original position:calculated for
%all cells
energymatrix.current_localPt=MF_SetNodePotential(1,LengthOfNode);
%figure;surfc(energymatrix.current_localPt)

nuclei_area=MF_GetNucleiFromVectData_img(CurrentModelMatrix);
CurrentModelMatrix.nuclei_area=nuclei_area;

%figure; imshow(nuclei_outlines,[]);impixelinfo
%figure; imshow(CurrentModelMatrix.nuclei_area,[]);impixelinfo


%%%%%%%%%%%%%evaluate nucleus / cell movement%%%%%%%%%%%%%%%%%%%%%
for iCell =ix
    if ~CurrentModelMatrix.FreezeTag(iCell) 

%the idea here is to find the preferential direction according to node
%position, some sort of leading edge thing.

movelefttoright=round(random('normal',0,1));
movetopdown=round(random('normal',0,1));

CurrentModelMatrix.Nuclei_Location(iCell,:)=CurrentModelMatrix.Nuclei_Location(iCell,:)+[movelefttoright,movetopdown];
energymatrix.cell_vector(iCell,:)=[movelefttoright,movetopdown];

%Evaluate whether the placing of the nuclei is done so that they are still
%enough separated.
%figure;imshow(CurrentModelMatrix.nuclei_area,[]);impixelinfo
CurrentModelMatrix.nuclei_area=MF_GetNucleiFromVectData_img(CurrentModelMatrix, iCell);
overlapping_list_Nuclei= MF_FindNucleiTension_list(CurrentModelMatrix,10,iCell);%rmk
%proceed with tension release as long as there are still overlapping nuclei

while ~isempty(overlapping_list_Nuclei) 
%determine which nuclei overlap
%release
CurrentModelMatrix=MF_ReleaseNucleiTensions(CurrentModelMatrix, overlapping_list_Nuclei);
listforcorrection=[overlapping_list_Nuclei(:,1);overlapping_list_Nuclei(:,2)]';
CurrentModelMatrix.nuclei_area=MF_GetNucleiFromVectData_img(CurrentModelMatrix, listforcorrection);
overlapping_list_Nuclei= MF_FindNucleiTension_list(CurrentModelMatrix,10,iCell);%rmk
end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Cells moved in module model_move_nodesSTD')
nuclei_outlines=MF_GetNucleiOutlinesFromVectData_img(CurrentModelMatrix);
disp('ok')

%until here everything is fine


for iCell =1:CurrentModelMatrix.numberofcells
    if ~CurrentModelMatrix.FreezeTag(iCell) 
%%%%Speed here could be improved
occupied_area=MF_GetCellAreaAroundSelectedNuclei_img(CurrentModelMatrix,interactionmatrix, iCell,true);


%calculate distribution of energy over nodes
%total_distribution=sum(energymatrix.distributionOfenergy{iCell}(:,1));

% PosKP = interactionmatrix.bordernodes{iCell,1}(:,5) == 0;
% DistancesFactorL = interactionmatrix.bordernodes{iCell,1}(PosKP,1);
% MinDistancesFactorL = interactionmatrix.bordernodes{iCell,1}(PosKP,6);
% 
% DiffDFL = DistancesFactorL - MinDistancesFactorL;
% [~,IX] = sort(DiffDFL); TMP = [IX interactionmatrix.DistribOfE]; Redistribfactors = sortrows(TMP,1);
% 
% energymatrix.distributionOfenergy{iCell}(:,1)= energymatrix.totalenergy(iCell,:).*Redistribfactors(:,2);%the longer the distance the more the energy.
% NORM = sum(energymatrix.distributionOfenergy{iCell}(:,1))./energymatrix.totalenergy(iCell,:);
% energymatrix.distributionOfenergy{iCell}(:,1) = energymatrix.distributionOfenergy{iCell}(:,1)./NORM;
% 
% 
% energymatrix.nodeenergy{iCell}(:,1)= energymatrix.distributionOfenergy{iCell}(:,1);%energymatrix.totalenergy(iCell,1)./total_distribution.*
cell_vector=energymatrix.cell_vector(iCell,:);
energymatrix.current_movementPt=MF_SetMovementPotential(cell_vector,1,LengthOfNode);


%figure; surfc(energymatrix.current_movementPt)
%potential as distance from nucleus and other cells: calculated per cell

%%%%this step must be improved%%%%
[energymatrix]=MF_SetWholeCellPotential(CurrentModelMatrix, interactionmatrix, iCell,occupied_area, ...
    CurrentModelMatrix.nuclei_area,nuclei_outlines,energymatrix,1,false);
%%%%
%figure;surfc((energymatrix.current_distancepotential+2*energymatrix.pt_other_cells+3*energymatrix.pt_nuclei_area))
%figure;surfc((energymatrix.current_distancepotential))
%figure;imshow(energymatrix.current_distancepotential,[]);impixelinfo
%figure;imshow((energymatrix.current_distancepotential+2*energymatrix.pt_other_cells+3*energymatrix.pt_nuclei_area),[]);impixelinfo

%%%%%%%%%%%%%%%%%%%do you really want to delete all nodes apart from the 7
%%%%%%%%%%%%%%%%%%%moving nodes?
interactionmatrix.bordernodes{iCell}(interactionmatrix.bordernodes{iCell}(:,5)>0,:)=[];
interactionmatrix.numberofbordernodes(iCell,1)=interactionmatrix.numberofnodes;
%nbNode=3
for nbNode=1:interactionmatrix.numberofnodes
% resulting_potential=zeros(31,31);

[resulting_potential, interactionmatrix]=MF_MergePotentials(CurrentModelMatrix,interactionmatrix,energymatrix, iCell,nbNode,LengthOfNode);

%figure;surf(resulting_potential);zlim([0 3])
%figure;imshow(resulting_potential,[]);impixelinfo
%choose random place below the value of the potential function at this
%place...
% potential_placesList=[];
potential_placesList = resulting_potential<energymatrix.nodeenergy{iCell}(nbNode,1);
%figure;imshow(resulting_potential,[]);impixelinfo
%testim=zeros(31,31);
%testim(potential_placesList)=1;
%figure;imshow(testim,[]);impixelinfo

if any(potential_placesList(:))
    potential_placesList = find(potential_placesList);
size_area_choosefrom=length(potential_placesList);


chosen_new_place_prob=round(random('Uniform',0.5,size_area_choosefrom+0.5));
chosen_new_place=potential_placesList(chosen_new_place_prob);
%testim2=zeros(31,31);
%testim2(chosen_new_place)=100;
%figure;imshow(testim2,[]);impixelinfo;
%resulting_potential(chosen_new_place);
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
end

coordinate=[];
coordinate(:,2)=interactionmatrix.bordernodes{iCell,1}(:,3);
coordinate(:,1)=interactionmatrix.bordernodes{iCell,1}(:,4);
%update now interactionmatrix for all nodes:
interactionmatrix.bordernodes{iCell,1}(:,1)=sqrt((CurrentModelMatrix.Nuclei_Location(iCell,1)-coordinate(:,2)).^2+ ...
    (CurrentModelMatrix.Nuclei_Location(iCell,2)-coordinate(:,1)).^2);
%calculate angle
angles=zeros(interactionmatrix.numberofnodes,1);
arctanindex=find(-CurrentModelMatrix.Nuclei_Location(iCell,1)+coordinate(:,2)>0);
arctanindexneg=find(-CurrentModelMatrix.Nuclei_Location(iCell,1)+coordinate(:,2)<0);
arctanindexzero=find(-CurrentModelMatrix.Nuclei_Location(iCell,1)+coordinate(:,2)==0);
angles(arctanindex)=atan((CurrentModelMatrix.Nuclei_Location(iCell,2)-coordinate(arctanindex,1))./ ...
    (-CurrentModelMatrix.Nuclei_Location(iCell,1)+coordinate(arctanindex,2)));
angles(arctanindexneg)=pi- atan((CurrentModelMatrix.Nuclei_Location(iCell,2)-coordinate(arctanindexneg,1))./ ...
    (CurrentModelMatrix.Nuclei_Location(iCell,1)-coordinate(arctanindexneg,2)));

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
end
%sort rows
interactionmatrix.bordernodes{iCell}=sortrows(interactionmatrix.bordernodes{iCell},2);

[~,interactionmatrix] = MF_BendMembrane(interactionmatrix,CurrentModelMatrix, iCell);
interactionmatrix.bordernodes{iCell}=sortrows(interactionmatrix.bordernodes{iCell},2);
    else
[interactionmatrix] = MF_ResetNodePosition(interactionmatrix,CurrentModelMatrix, iCell,1);        
    end
end %over all cells

disp('New Node position evaluated in module model_move_nodesSTD')

%virtualimage=visualize_celloutline_model_cell(interactionmatrix,CurrentModelMatrix);
%virtualimage(virtualimage==8)=200;
%figure;imshow(virtualimage,[]);impixelinfo
%temp1=virtualimage(1100:1500,1600:2150);
%actualimage=figure;imshow(temp1,[])
%filename=sprintf('C:\\Users\\Stephan Daetwyler\\Documents\\ETH\\BMasterwork\\Simulations\\movement2\\actualimage%d.png',iteration_nb)
%print(actualimage, '-dpng', filename);
