function[readytodivide,CurrentModelMatrix]=MF_WhichCellMustDivide(CurrentModelMatrix)%,virtualimage)
%[MF]Massive changes 06/13
%This module gives a tag for each cell, a random value in hour between 20 and 24.
%when the lifetime of a cell reaches its tag, the cell is ready to divide



%%%%New part to randomize only once, which is more realistic behaviour%%%
if ~isfield(CurrentModelMatrix,'WhenToDivide')
    CurrentModelMatrix.WhenToDivide = normrnd(22,2,CurrentModelMatrix.numberofcells,1);
elseif size(CurrentModelMatrix.WhenToDivide,1) < CurrentModelMatrix.numberofcells
    extraCells = normrnd(22,2,(CurrentModelMatrix.numberofcells - size(CurrentModelMatrix.WhenToDivide,1)),1);
    CurrentModelMatrix.WhenToDivide = [CurrentModelMatrix.WhenToDivide;extraCells];
end

readytodivide=find(CurrentModelMatrix.status>CurrentModelMatrix.WhenToDivide & ~CurrentModelMatrix.FreezeTag).';%turns LCD in exponential, why?%+ exp(ModelMatrix.LCD) normrnd(22,2,ModelMatrix.numberofcells,1)


 
%show which cells were selected for division
%colorimage1 =zeros(ModelMatrix.rownumber,ModelMatrix.columnnumber) ;
%colorimage2 =zeros(ModelMatrix.rownumber,ModelMatrix.columnnumber) ;
%colorimage3 =zeros(ModelMatrix.rownumber,ModelMatrix.columnnumber) ;

%for i = [isx]
 %   ix = find(virtualimage ==ModelMatrix.ObjectID(i));
  %  colorimage1(ix)=ModelMatrix.color(i,1);
   % colorimage2(ix)=ModelMatrix.color(i,2);
    %colorimage3(ix)=ModelMatrix.color(i,3);
%end

%vir=ceil(mat2gray(virtualimage))*0.5+colorimage1;
%colorimage=cat(3,vir,colorimage2,colorimage3);
%imshow(colorimage)  

