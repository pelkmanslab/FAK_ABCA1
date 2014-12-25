function[fast_outline,nuclei_area]=MF_GetNucleiOutlinesFromVectData_img(CurrentModelMatrix)
%[MF] check and optimized 06/13, not much
%%%%%%%%%%%%%%%%%%%get nuclei outlines
nuclei_area=MF_GetNucleiFromVectData_img(CurrentModelMatrix);
fast_outline = nuclei_area;
Pos = bwmorph(nuclei_area,'remove');
Pos = ~Pos;
fast_outline(Pos) = 0;






