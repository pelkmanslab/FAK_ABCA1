function [interactionmatrix,CurrentModelMatrix,energymatrix,ReRun] = MF_SyncMovementOfNodes(energymatrix,CurrentModelMatrix,interactionmatrix,nodeGrowthlimit,redo)

ReRun = false;
[overlaplist, interactionmatrix]=MF_FindOverlapsAndGuiltyNodes(CurrentModelMatrix, interactionmatrix);
% figure(2);imshow(IMGwithOverlap)
%second release
count = 0;
while ~isempty(overlaplist)
    count = count+1;%test = (1:1000)';
%      if any((count./5) == test)
%          param = [0,1];
% %      elseif any((count./10)==test)
%      else
         param = [0,0];
%      end
     
     
     if  count == redo
        ReRun = true;
        %disp('redo cell division step')
        break
     end
    
    tic
    [interactionmatrix,CurrentModelMatrix,energymatrix]=...
        MF_SyncNodesOnly(param,energymatrix,CurrentModelMatrix,interactionmatrix, overlaplist,nodeGrowthlimit);
    [overlaplist, interactionmatrix]=MF_FindOverlapsAndGuiltyNodes(CurrentModelMatrix, interactionmatrix);
    %figure(1);imshow(IMGwithOverlap)
    toc 
    
    
    if ~isempty(overlaplist)%any((count./10) == test)
    [CurrentModelMatrix, interactionmatrix,energymatrix]=MF_TOcheck(CurrentModelMatrix, interactionmatrix,energymatrix);  
    [overlaplist, interactionmatrix]=MF_FindOverlapsAndGuiltyNodes(CurrentModelMatrix, interactionmatrix);
    %figure(1);imshow(IMGwithOverlap)
    end
    
    
    if ~isempty(overlaplist)
    [interactionmatrix,CurrentModelMatrix,energymatrix]=MF_SyncNodesWithRandNuclWiggling(energymatrix,CurrentModelMatrix,interactionmatrix,nodeGrowthlimit);
    [overlaplist, interactionmatrix]=MF_FindOverlapsAndGuiltyNodes(CurrentModelMatrix, interactionmatrix);
    end
end
end