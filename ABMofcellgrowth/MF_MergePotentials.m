function[sum_of_potentials,interactionmatrix]=MF_MergePotentials(CurrentModelMatrix,interactionmatrix,energymatrix, iCell,nbNode,LengthOfNode)
%[MF] Not satisfying 06/13
%sum of potential for a single node


% smoothNucleiEnergy = CurrentModelMatrix.smoothNucleiEnergy;
% if size(smoothNucleiEnergy,1) ~= size(energymatrix.pt_nuclei_area,1) || ...
%       size(smoothNucleiEnergy,2) ~= size(energymatrix.pt_nuclei_area,2) 
%   logicalIMG = energymatrix.pt_nuclei_area == max(energymatrix.pt_nuclei_area(:));
%   TMP = regionprops(logicalIMG,'centroid');TMP =round(TMP.Centroid);
%   Shift = TMP - [125 125];
%   if Shift(1,1)<0
%   InterrestingX = [abs(Shift(1,1)) 251]; AddTMPAfterX = size(energymatrix.pt_nuclei_area,1) - (251 - abs(Shift(1,1)) + 1);
%   else
%   InterrestingX = [1 251-abs(Shift(1,1))]; AddTMPBeforeX = size(energymatrix.pt_nuclei_area,1) - (251 - abs(Shift(1,1)));
%   end
%   if Shift(1,2)<0
%   InterrestingY = [abs(Shift(1,2)) 251]; AddTMPAfterY = size(energymatrix.pt_nuclei_area,2) - (251 - abs(Shift(1,2)) + 1);
%   else
%   InterrestingY = [1 251-abs(Shift(1,2))]; AddTMPBeforeY = size(energymatrix.pt_nuclei_area,2) - (251 - abs(Shift(1,2)));
%   end
%   TMP = smoothNucleiEnergy(InterrestingX(1,1):InterrestingX(1,2),InterrestingY(1,1):InterrestingY(1,2));
%   
%   if Shift(1,1)<0 && Shift(1,2)<0
%   TMP = [TMP;zeros(AddTMPAfterX,size(TMP,2))];
%   TMP = [TMP zeros(size(TMP,1),AddTMPAfterY)];
%   elseif Shift(1,1)>0 && Shift(1,2)>0
%   TMP = [zeros(AddTMPBeforeX,size(TMP,2));TMP];
%   TMP = [zeros(size(TMP,1),AddTMPBeforeY) TMP];  
%   elseif Shift(1,1)<0 && Shift(1,2)>0
%   TMP = [TMP;zeros(AddTMPAfterX,size(TMP,2))];
%   TMP = [zeros(size(TMP,1),AddTMPBeforeY) TMP];    
%   elseif Shift(1,1)>0 && Shift(1,2)<0
%   TMP = [zeros(AddTMPBeforeX,size(TMP,2));TMP];
%   TMP = [TMP zeros(size(TMP,1),AddTMPAfterY)];
%   elseif Shift(1,1)<0 && Shift(1,2)== 0
%   TMP = [TMP;zeros(AddTMPAfterX,size(TMP,2))];
%   elseif Shift(1,1)>0 && Shift(1,2)== 0
%   TMP = [zeros(AddTMPBeforeX,size(TMP,2));TMP];    
%   elseif Shift(1,1)==0 && Shift(1,2)< 0
%   TMP = [TMP zeros(size(TMP,1),AddTMPAfterY)];
%   elseif Shift(1,1)==0 && Shift(1,2)> 0 
%   TMP = [zeros(size(TMP,1),AddTMPBeforeY) TMP];
%   end
% smoothNucleiEnergy = TMP;  
% clear TMP
% end
% NormNenergy = energymatrix.pt_nuclei_area./max(energymatrix.pt_nuclei_area(:));
% kernel = ones(15, 15) ; % 3x3 mean kernel
smoothNucleiEnergy = energymatrix.pt_nuclei_area;
% smoothNucleiEnergy = conv2(smoothNucleiEnergy, kernel, 'same'); % Convolve keeping size of Img
smoothNucleiEnergy = smoothNucleiEnergy./max(smoothNucleiEnergy(:));
% smoothNucleiEnergy(NormNenergy == 1) = 1;
wholecellpotential= energymatrix.current_distancepotential+4*logical(energymatrix.pt_other_cells)+smoothNucleiEnergy; %made such for nodes to not be able to cross nuclei or crawl over another cell
%figure;imshow(temp2,[]);impixelinfo
movementPt=energymatrix.current_movementPt./10;
localPt=energymatrix.current_localPt./10;
%temp=energymatrix.currentnodelocations;
%temp(temp>0)=0.8;
%temp2=temp+wholecellpotential;
%imshow(temp,[])
%figure;imshow(energymatrix.currentnodelocations,[]);impixelinfo
%figure;mesh(wholecellpotential);impixelinfo

%idea find pixels in image with pixel nbs to perfectly align overlay of
%potential functions around pixels.
%Deal also with the problem that two or more pixels might be placed at the same place. 

cond1=interactionmatrix.bordernodes{iCell,1}(nbNode,4)<=1 ...
    | interactionmatrix.bordernodes{iCell,1}(nbNode,4)>=CurrentModelMatrix.rownumber ...
    | interactionmatrix.bordernodes{iCell,1}(nbNode,3)<=1 ...
    | interactionmatrix.bordernodes{iCell,1}(nbNode,3)>=CurrentModelMatrix.columnnumber; 

if ~isempty(find(energymatrix.currentnodelocations==10000+nbNode,1))
    [temppos,temppos2]=find(energymatrix.currentnodelocations==10000+nbNode);
 elseif cond1 %if it is not recognized e.g. at the border or other strange bugs...
         index_temp1=find(energymatrix.currentnodelocations>10000);
         
      
         index_temp1=index_temp1(1);
         value_known_node=energymatrix.currentnodelocations(index_temp1)-10000;
         [temppos,temppos2]=find(energymatrix.currentnodelocations==10000+value_known_node);
     

         %figure;imshow(energymatrix.currentnodelocations,[]);impixelinfo
elseif interactionmatrix.bordernodes{iCell,1}(nbNode,1)>3 * interactionmatrix.bordernodes{iCell,1}(nbNode,6)
    
    interactionmatrix.bordernodes{iCell,1}(nbNode,1)=interactionmatrix.bordernodes{iCell,1}(nbNode,6)+3;
    dist1=interactionmatrix.bordernodes{iCell,1}(nbNode,1)*cos(interactionmatrix.bordernodes{iCell,1}(nbNode,2));
    dist2=interactionmatrix.bordernodes{iCell,1}(nbNode,1)*sin(interactionmatrix.bordernodes{iCell,1}(nbNode,2));
    interactionmatrix.bordernodes{iCell,1}(nbNode,3)=round(CurrentModelMatrix.Nuclei_Location(iCell,1)+dist1);
    interactionmatrix.bordernodes{iCell,1}(nbNode,4)=round(CurrentModelMatrix.Nuclei_Location(iCell,2)-dist2);

         
         index_temp1=find(energymatrix.currentnodelocations>10000);
     
         index_temp1=index_temp1(1);
         value_known_node=energymatrix.currentnodelocations(index_temp1)-10000;
         [temppos,temppos2]=find(energymatrix.currentnodelocations==10000+value_known_node);

 else %deal with case that two or more nodes are placed at the same pixel!
        istem2=[];
        for ind=1:interactionmatrix.numberofbordernodes(iCell)
        if ~isempty(find(energymatrix.currentnodelocations==10000+ind,1))
        istem2=[istem2,ind];
        end
        end
        
        indextemp=min(istem2(istem2>nbNode));
        if ~isempty(indextemp)
       [temppos,temppos2]=find(energymatrix.currentnodelocations==10000+indextemp);
        else
        indextemp=max(istem2);  
       [temppos,temppos2]=find(energymatrix.currentnodelocations==10000+indextemp);
        end
end

relPos=[temppos2,temppos];

clear temppos2 temppos indextemp i istem2 istem1 ind

%figure;imshow(energymatrix.currentnodelocations,[]);impixelinfo

if ~isempty(relPos)
small_window=LengthOfNode;
try%%%tatata this is sooo bad.
zoom_in_node=wholecellpotential((relPos(1,2)-small_window):(relPos(1,2)+small_window),(relPos(1,1)-small_window):(relPos(1,1)+small_window));
catch
    [temp1,temp2]=size(wholecellpotential);
    addedpotential=zeros(temp1+80,temp2+80);
    addedpotential(:,:)=9;
    addedpotential(40:(temp1+39),(40:temp2+39))=wholecellpotential;
    relPos(1,1)=relPos(1,1)+39;
    relPos(1,2)=relPos(1,2)+39;
    zoom_in_node=addedpotential((relPos(1,2)-small_window):(relPos(1,2)+small_window),(relPos(1,1)-small_window):(relPos(1,1)+small_window));
end
else
    zoom_in_node = 0;
end 
    
%figure;imshow(zoom_in_node,[]);impixelinfo
%energymatrix.factors=[1,1,1];
sum_of_potentials=energymatrix.factors(1,1)*zoom_in_node+energymatrix.factors(1,2)*localPt+energymatrix.factors(1,3)*movementPt;

%figure;imshow(sum_of_potentials,[]);impixelinfo
%figure;surf(movementPt)

%plot where nodes will be placed relative to center...
%figure;imshow(zoom_in_node,[]);impixelinfo
%outpict=zeros(201,201);
%relind=relPos1(:,1).*201+relPos1(:,2);
%outpict(relind)=1;
%figure;imshow(outpict,[]);impixelinfo
