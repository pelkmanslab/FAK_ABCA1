function[distancePotential]=MF_SetNodePotential(normalized,LengthOfNode)
%draw potential for distance away from old point:calculated for all cells
%in the same way. If normalized=1, then normalize output function

%test:normalized=1

BW=zeros(31,31);
BW(16,16)=1;
distancefromnode = bwdist(BW);
distancePotential=double(distancefromnode.^2);

if LengthOfNode <= 15
distancePotential = distancePotential(16-LengthOfNode:16+LengthOfNode,16-LengthOfNode:16+LengthOfNode);
end

if normalized == 1
   distancePotential=distancePotential/max(distancePotential(:));
end
   
