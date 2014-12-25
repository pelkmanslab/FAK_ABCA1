function [CurrentModelMatrix, interactionmatrix] = MF_ReleaseDivCellTension(CurrentModelMatrix,interactionmatrix, overlapping_list)
%[MF] massive changes on 06/13



if ~isempty(overlapping_list)
    for index = 1: length(overlapping_list(:,1))
    nuc1=overlapping_list(index,1);
    nuc2=overlapping_list(index,2);
    ang=overlapping_list(index,3);
  
% %   if CurrentModelMatrix.growthtimestepincrement < 1
%     distr1 = random('normal',14,3)+(10*CurrentModelMatrix.growthtimestepincrement);%be careful this value has to be always more than the threshold for saying that there is an overlap
%     distr2= random('normal',14,3)+(10*CurrentModelMatrix.growthtimestepincrement); 
%   else
    distr1 =random('normal',3,0.5);
    distr2= random('normal',3,0.5);
%   end
    shift=[distr1*cos(ang),-distr1*sin(ang)];
    negshift=[distr2*cos(ang+pi),-distr2*sin(ang+pi)];
    
    %be careful that the tension release does not exceeds matrix dimensions
    cond1=CurrentModelMatrix.Nuclei_Location(nuc1,1)+negshift(1);
    cond2=CurrentModelMatrix.Nuclei_Location(nuc1,2)+negshift(2);
    cond3=CurrentModelMatrix.Nuclei_Location(nuc2,1)+shift(1);
    cond4=CurrentModelMatrix.Nuclei_Location(nuc2,2)+shift(2);
    
   %some problems can appear here, let's try to make it smarter. 
   %as it is now, when it comes to the border, it just stop moving, so be
   %it, but then it must not divide anymore, simulating divisions occuring in the
   %unseen area
   if cond1>0 && cond2>0 && cond1<CurrentModelMatrix.columnnumber && cond2<CurrentModelMatrix.rownumber && CurrentModelMatrix.FreezeTag(nuc1,1) == false
    CurrentModelMatrix.Nuclei_Location(nuc1,:)=CurrentModelMatrix.Nuclei_Location(nuc1,:) + negshift;
   else
    CurrentModelMatrix.FreezeTag(nuc1,1) = true;  
   end
   if cond3>0 && cond4>0 && cond3<CurrentModelMatrix.columnnumber && cond4<CurrentModelMatrix.rownumber && CurrentModelMatrix.FreezeTag(nuc2,1) == false
    CurrentModelMatrix.Nuclei_Location(nuc2,:) =CurrentModelMatrix.Nuclei_Location(nuc2,:) + shift;
   else
    CurrentModelMatrix.FreezeTag(nuc2,1) = true;    
   end
   %now cells that must not move and shouldnt get touched anymore are marked with a 1
   %in the .FreezeTag list(logical)
   
   %move the nodes accordingly
   if CurrentModelMatrix.FreezeTag(nuc1,1) == false;
       NodeData = interactionmatrix.bordernodes{nuc1};
       NodeData(:,3)= round(NodeData(:,3)+negshift(1));
       NodeData(:,4)= round(NodeData(:,4)+negshift(2));
       %update
       NodeData(:,1)=sqrt((CurrentModelMatrix.Nuclei_Location(nuc1,1)-NodeData(:,3)).^2+ ...
       (CurrentModelMatrix.Nuclei_Location(nuc1,2)-NodeData(:,4)).^2);
       interactionmatrix.bordernodes{nuc1} = NodeData;
   end
   
   if CurrentModelMatrix.FreezeTag(nuc2,1) == false;
       NodeData = interactionmatrix.bordernodes{nuc2};
       NodeData(:,3)= round(NodeData(:,3)+shift(1));
       NodeData(:,4)= round(NodeData(:,4)+shift(2));
       %update
       NodeData(:,1)=sqrt((CurrentModelMatrix.Nuclei_Location(nuc2,1)-NodeData(:,3)).^2+ ...
       (CurrentModelMatrix.Nuclei_Location(nuc2,2)-NodeData(:,4)).^2);
       interactionmatrix.bordernodes{nuc2} = NodeData;
   end
   
    end
    
end

