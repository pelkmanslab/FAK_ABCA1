function [CurrentModelMatrix,interactionmatrix,energymatrix] = MF_MakeSpace4DividingCells(CurrentModelMatrix,interactionmatrix,energymatrix,isx,trshld)

% LOCTGbis = MF_FindNucleiTension_list(CurrentModelMatrix,trshld,isx);
% while ~isempty(LOCTGbis)
listOfcelltoGivespace = MF_FindNucleiTension_list(CurrentModelMatrix,trshld,isx);
while ~isempty(listOfcelltoGivespace)
%first make space around the nuclei that is about to divide.
[CurrentModelMatrix,interactionmatrix]=MF_ReleaseDivCellTension(CurrentModelMatrix,interactionmatrix, listOfcelltoGivespace);
%determine which nuclei still overlap

%node overlap fixing

[interactionmatrix,CurrentModelMatrix,energymatrix] = MF_SyncMovementOfNodes(energymatrix,CurrentModelMatrix,interactionmatrix,60,inf);

listOfcelltoGivespace = MF_FindNucleiTension_list(CurrentModelMatrix,trshld,isx);

end
%then release normal tensions
% overlapping_list= MF_FindNucleiTension_list(CurrentModelMatrix,10,isx);
% while ~isempty(overlapping_list)
% %release
% [CurrentModelMatrix,interactionmatrix]=MF_ReleaseDivCellTension(CurrentModelMatrix,interactionmatrix, overlapping_list);
% %determine which nuclei overlap
% overlapping_list= MF_FindNucleiTension_list(CurrentModelMatrix,8,isx);
% 
% [interactionmatrix,CurrentModelMatrix,energymatrix] = MF_SyncMovementOfNodes(energymatrix,CurrentModelMatrix,interactionmatrix,60,inf);
% 
% end
% %then recheck for the specific nuclei that will divide, if it is not good,
% %redo!
% LOCTGbis = MF_FindNucleiTension_list(CurrentModelMatrix,trshld,isx);
% end
end