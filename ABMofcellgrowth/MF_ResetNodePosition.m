function [interactionmatrix] = MF_ResetNodePosition(interactionmatrix,CurrentModelMatrix, iCell,offset)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%define a set of new nodes.
cutStep = interactionmatrix.numberofnodes;
incrementstep=2*pi/cutStep;
angleinnernode = 0:incrementstep:(2*pi - incrementstep); %hardcoded!
howmanyinnernodes=length(angleinnernode);
interactionmatrix.bordernodes{iCell,1}=zeros(interactionmatrix.numberofnodes,6);
interactionmatrix.numberofbordernodes(iCell)=interactionmatrix.numberofnodes;
%calculate here how large is the distance to nuclei membrane from centroid  
factorl=zeros(howmanyinnernodes,1);
innernodeXpos=zeros(howmanyinnernodes,1);
innernodeYpos=zeros(howmanyinnernodes,1);

%defines positions of the innernodes, reminder for those who needs it: if
%A=tan(B) then B=atan(A). You are welcome ;)
for innernode=1:howmanyinnernodes
    factort=atan(CurrentModelMatrix.radius(iCell,1)/CurrentModelMatrix.radius(iCell,2) ...
        *tan(angleinnernode(innernode)-CurrentModelMatrix.angle(iCell,1)));
    factorl(innernode,1)=sqrt((CurrentModelMatrix.radius(iCell,1)*cos(factort))^2 + (CurrentModelMatrix.radius(iCell,2)*sin(factort))^2);
    factorl(innernode,1)=factorl(innernode,1)+offset;%hardcoded too !!!
    dist1=factorl(innernode)*cos(angleinnernode(innernode));
    dist2=factorl(innernode)*sin(angleinnernode(innernode));
innernodeXpos(innernode,1)=round(CurrentModelMatrix.Nuclei_Location(iCell,1)+dist1);
innernodeYpos(innernode,1)=round(CurrentModelMatrix.Nuclei_Location(iCell,2)-dist2);
end   
%we now have all coordinates and angles and norms of the new nodes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   interactionmatrix.bordernodes{iCell}(:,1)=factorl;
   interactionmatrix.bordernodes{iCell}(:,2)=angleinnernode';
   interactionmatrix.bordernodes{iCell}(:,3)=innernodeXpos;
   interactionmatrix.bordernodes{iCell}(:,4)=innernodeYpos;
   interactionmatrix.bordernodes{iCell}(:,5)=0;
   interactionmatrix.bordernodes{iCell}(:,6)=factorl - offset + 1;
   interactionmatrix.bordernodes{iCell}=sortrows(interactionmatrix.bordernodes{iCell},2);
%update now angle_nodesofcell for calculation (several innernodes might be
%placed per cell in this loop 
   
   
end

 
 %visualize here the outlines of the cell
 %virtualimage=visualize_celloutline_model_cell(interactionmatrix,CurrentModelMatrix);
   
   

