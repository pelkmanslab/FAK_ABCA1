function[distancefromolddirection_potential]=MF_SetMovementPotential(cell_vector,normalized,LengthOfNode)
%[MF] fine 06/13 need to comment more
%bias towards whole cell movement (cell_vector):calculated per cell
%test: cell_vector=[1,0]

special_case_00=0;
%move top down
if cell_vector(1,2)~=0
    if cell_vector(1,2)<0
        lambda=-16/cell_vector(1,2);
        pos1=1;
        pos2=16+lambda*cell_vector(1,1);%move right left
    else
        lambda=16/cell_vector(1,2);
        pos1=30;
        pos2=15+lambda*cell_vector(1,1);
    end
else
    if cell_vector(1,1)>0
        pos1=16;
        pos2=30;
    elseif cell_vector(1,1)<0
        pos1=16;
        pos2=0;
    elseif cell_vector(1,1)==0
        special_case_00=1;
    end
end

if special_case_00==1
    temp1=zeros(31,31);
    temp1(16,16)=1;
else
temp1=linept(zeros(31,31),16,16,pos1,pos2);
temp1=temp1(1:31,1:31);
temp1(temp1>0)=1;
end

distancefromolddirection=double(bwdist(temp1));
distancefromolddirection_potential=(distancefromolddirection.^2);

if LengthOfNode <= 15
distancefromolddirection_potential = distancefromolddirection_potential(16-LengthOfNode:16+LengthOfNode,16-LengthOfNode:16+LengthOfNode);
end

if normalized==1
distancefromolddirection_potential=distancefromolddirection_potential/max(distancefromolddirection_potential(:));
end

