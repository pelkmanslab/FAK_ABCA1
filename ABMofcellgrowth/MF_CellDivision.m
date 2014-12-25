function[CurrentModelMatrix,interactionmatrix,energymatrix]=MF_CellDivision(CurrentModelMatrix,interactionmatrix,energymatrix,virtualimage,isx,totalenergy, minCellsize)
%[MF] minor changes 06/13, principally the skipping of  frozen cells, a bit
%of polishing, and heavy commenting
%Cell division Control panel, frozen cells are not processed, they are
%filtered out in the whatstagecellis module, updated in the release tension
%function

energymatrix.HasitJustDivided = false(CurrentModelMatrix.numberofcells,1);

%1. Check which cells should divide:
if ~isempty(isx)
    
%make space for the cell to divide
%the result of these precise nested while loops is that there will be
%enough space just around the dividing nuclei, this will spread through
%populations of cells, by a neighbour to neighbour effect that will make
%islets growing ultimatelly. this is mimicking exactly what happens in a
%culture of cells

      

listOfCells = [];
for allcells = isx %there might be more than one cell dividing per frame (loop over all those)
numberofcells=CurrentModelMatrix.numberofcells;


%Update all information in CurrentModelMatrix (Nuclei Information)
[CurrentModelMatrix,energymatrix,~]=MF_CellDivisionNodePlacement(CurrentModelMatrix, virtualimage, allcells, minCellsize,interactionmatrix,energymatrix);%ok checked

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Move new Nuclei apart
overlapping_list= MF_FindNucleiTension_list(CurrentModelMatrix);
%proceed with tension release as long as there are still overlapping nuclei
 while ~isempty(overlapping_list)
%release
CurrentModelMatrix=MF_ReleaseNucleiTensions(CurrentModelMatrix, overlapping_list);
%determine which nuclei overlap
overlapping_list= MF_FindNucleiTension_list(CurrentModelMatrix);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%place now new nodes around the cell
interactionmatrix.bordernodes{allcells}=[];%reset all nodes
interactionmatrix.bordernodes{numberofcells+1}=[];

for whichcell=[allcells, numberofcells+1]   
    
    listOfCells = [listOfCells whichcell];
%define a set of new nodes. 
angle=random('Uniform',0,(2*pi),1,interactionmatrix.numberofnodes);
angle=sort(angle);

%calculate here how large is the distance centroid-nuclei membrane   
factorl=zeros(length(angle),1);
nodeXpos=zeros(length(angle),1);
nodeYpos=zeros(length(angle),1);

size1 = 80*CurrentModelMatrix.growthtimestepincrement;
size2 = 6*CurrentModelMatrix.growthtimestepincrement;

%set coordinates
for newnode=1:length(angle)    
    factort=atan(CurrentModelMatrix.radius(whichcell,1)/CurrentModelMatrix.radius(whichcell,2)*tan(angle(newnode)-CurrentModelMatrix.angle(whichcell,1)));
    factorl(newnode,1)=sqrt((CurrentModelMatrix.radius(whichcell,1)*cos(factort))^2 + (CurrentModelMatrix.radius(whichcell,2)*sin(factort))^2);
    factorl(newnode,1)=factorl(newnode,1) + random('normal',size1,size2);%randomness for the new position, Mustbefunction of the timescale
    dist1=factorl(newnode)*cos(angle(newnode));
    dist2=factorl(newnode)*sin(angle(newnode));
nodeXpos(newnode,1)=round(CurrentModelMatrix.Nuclei_Location(whichcell,1)+dist1);
nodeYpos(newnode,1)=round(CurrentModelMatrix.Nuclei_Location(whichcell,2)+dist2);
end   

interactionmatrix.bordernodes{whichcell,1}(:,1)=factorl;%norm
interactionmatrix.bordernodes{whichcell,1}(:,2)=angle;
interactionmatrix.bordernodes{whichcell,1}(:,3)=nodeXpos;
interactionmatrix.bordernodes{whichcell,1}(:,4)=nodeYpos;
interactionmatrix.bordernodes{whichcell,1}(:,5)=0;
interactionmatrix.bordernodes{whichcell,1}(:,6)=factorl;%.... let see why he did this, but erase is pending

interactionmatrix.numberofbordernodes(whichcell,1)=interactionmatrix.numberofnodes;
 

[~,interactionmatrix] = MF_BendMembrane(interactionmatrix,CurrentModelMatrix, whichcell);%frozen cells are already not taken into account in the loop 


%here it becomes interresting ;) we set the energy environment
energymatrix.cell_vector(whichcell,:)=[0,0];
energymatrix.totalenergy(whichcell,:)=random('normal',totalenergy(1,1),totalenergy(1,2));

%lets give more energy to nodes that are far and less to nodes that are
%close
PosKP = interactionmatrix.bordernodes{whichcell,1}(:,5) == 0;
DistancesFactorL = interactionmatrix.bordernodes{whichcell,1}(PosKP,1);
MinDistancesFactorL = interactionmatrix.bordernodes{whichcell,1}(PosKP,6);
DiffDFL = DistancesFactorL - MinDistancesFactorL;
[~,IX] = sort(DiffDFL); TMP = [IX interactionmatrix.DistribOfE]; Redistribfactors = sortrows(TMP,1);

energymatrix.distributionOfenergy{whichcell}(:,1) = (energymatrix.totalenergy(whichcell,:)./interactionmatrix.numberofnodes).*Redistribfactors(:,2);%the longer the distance the more the energy.
NORM = sum(energymatrix.distributionOfenergy{whichcell}(:,1))./energymatrix.totalenergy(whichcell,:);
energymatrix.distributionOfenergy{whichcell}(:,1) = energymatrix.distributionOfenergy{whichcell}(:,1)./NORM;

%total_distribution=sum(energymatrix.distributionOfenergy{whichcell}(:,1));
energymatrix.nodeenergy{whichcell}(:,1)= energymatrix.distributionOfenergy{whichcell}(:,1);%energymatrix.totalenergy(whichcell,1)./total_distribution.*

energymatrix.HasitJustDivided(whichcell) = true;


end
 CurrentModelMatrix.numberofcells = numberofcells + 1;
 

end




[nuclei_outlines,~]=MF_GetNucleiOutlinesFromVectData_img(CurrentModelMatrix);
occupied_area=MF_GetCellArea_img(CurrentModelMatrix,interactionmatrix,true);
CurrentModelMatrix.nuclei_area=MF_GetNucleiFromVectData_img(CurrentModelMatrix);
%imshow(occupied_area )
 for whichcell = listOfCells
    %define the accessible area for the new cells now that the objects
    %exists, I do it now as all structures have been updated. It avoids anoying bugs.   
 
 
 energymatrix = MF_SetWholeCellPotential(CurrentModelMatrix, interactionmatrix, whichcell,occupied_area, ...
 CurrentModelMatrix.nuclei_area,nuclei_outlines,energymatrix,1,false);

 accessibleArea = energymatrix.pt_allsurroundingNucl + logical(energymatrix.pt_other_cells) + (2.*energymatrix.pt_nuclei_area); %this allows me to put the nodes at the right place.
 %imshow(accessibleArea)
 %first replace the nodes in the accessible area.
 PosN2BC = 1;
 PosNinNucl = 1;
 nodesCaract = interactionmatrix.bordernodes{whichcell};
 nodeX = nodesCaract(:,4);replacedX =floor( nodeX - CurrentModelMatrix.Nuclei_Location(whichcell,2) + CurrentModelMatrix.EnRadius + 1);
 nodeY = nodesCaract(:,3);replacedY =floor( nodeY - CurrentModelMatrix.Nuclei_Location(whichcell,1) + CurrentModelMatrix.EnRadius + 1);
 
 %if there is any node placed
 
 
 
 
 OriginalNorm = nodesCaract(:,1);
 disp('control of dividing cell placement, careful')
 count = 0;
 while any(PosN2BC) || any(PosNinNucl)   
 count = count + 1;
 if count > 10000 
     break
 end
 OUT = replacedY < 1 | replacedY > ((CurrentModelMatrix.EnRadius.*2) + 1) | replacedX < 1 | replacedX > ((CurrentModelMatrix.EnRadius.*2) + 1);
 absolutePos = ((CurrentModelMatrix.EnRadius.*2) + 1).*(replacedY-1) + replacedX;%which explains the inversion of x and Y here.
 absolutePos(OUT) = 1; %safety, indeed nodes can be out from the area, as we bring them back, they will reenter in the area.
 PosN2BC = accessibleArea(absolutePos) == 1;%gets a one if in another cell.
 PosN2BC(OUT)= true;
 PosNinNucl = accessibleArea(absolutePos) == 2 & nodesCaract(:,1) <= (OriginalNorm+1); %gets a one if in its nuclei and if the norm has not reached the initial size
 nodesToBringCloser = nodesCaract(PosN2BC,1);
 nodesToPutFarther = nodesCaract(PosNinNucl,1);
 nodesToBringCloser = nodesToBringCloser - 1 ;
 nodesToPutFarther = nodesToPutFarther + 1 ; 
 %let's update the coordinates
 nodesCaract(PosN2BC,1) = nodesToBringCloser;
 nodesCaract(PosNinNucl,1) = nodesToPutFarther;
 replacedX = round(nodesCaract(:,1).*sin(nodesCaract(:,2)) + CurrentModelMatrix.EnRadius + 1);
 replacedY = round(nodesCaract(:,1).*cos(nodesCaract(:,2)) + CurrentModelMatrix.EnRadius + 1);
 %think about a safety escape
 end
 %nodes are placed correctly, let's replace the node in the real image
 
 nodesCaract(:,4) = round(replacedX - (CurrentModelMatrix.EnRadius + 1) + CurrentModelMatrix.Nuclei_Location(whichcell,2)); 
 nodesCaract(:,3) = round(replacedY - (CurrentModelMatrix.EnRadius + 1) + CurrentModelMatrix.Nuclei_Location(whichcell,1));
 disp('everything went fine')
 interactionmatrix.bordernodes{whichcell} = nodesCaract;
 
 
[nuclei_outlines,~]=MF_GetNucleiOutlinesFromVectData_img(CurrentModelMatrix);
occupied_area=MF_GetCellArea_img(CurrentModelMatrix,interactionmatrix,true);
CurrentModelMatrix.nuclei_area=MF_GetNucleiFromVectData_img(CurrentModelMatrix);
 
 end
end

