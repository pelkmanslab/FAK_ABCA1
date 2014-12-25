function[virtualimage]=MF_GetNucleiFromVectData_img(CurrentModelMatrix,iCell)



if nargin == 2
for i=iCell
    
CurrentModelMatrix.nuclei_area(CurrentModelMatrix.nuclei_area==CurrentModelMatrix.ObjectID(i))=0;
    
end
virtualimage = CurrentModelMatrix.nuclei_area;
totalcellnumber = iCell;
else
totalcellnumber = 1:CurrentModelMatrix.numberofcells;
virtualimage = uint16(zeros(CurrentModelMatrix.rownumber,CurrentModelMatrix.columnnumber));
end

%%%%%%%%%%%%%Calculations for Nucleus%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
center = round(CurrentModelMatrix.Nuclei_Location);
radius = round(CurrentModelMatrix.radius);
columnnumber = CurrentModelMatrix.columnnumber;
rownumber = CurrentModelMatrix.rownumber;

%defines border values of image which is cut out 
leftboundary  = center(:,1)-radius(:,1);
leftboundary(leftboundary<1)=1;
rightboundary = center(:,1)+radius(:,1);
rightboundary(rightboundary>columnnumber)=columnnumber;
upperboundary = center(:,2)-radius(:,2);
upperboundary(upperboundary<1)=1;
lowerboundary = center(:,2)+radius(:,2);
lowerboundary(lowerboundary>rownumber)=rownumber;


for i = totalcellnumber
   
%initializes rotation
theta=CurrentModelMatrix.angle(i);
R=[cos(theta), -sin(theta),0; sin(theta), cos(theta),0;0,0,1];
tform=maketform('affine',R);

%makes an ellipse

X = ones(2*radius(i,2)+1,1)*(-radius(i,1):(radius(i,1)));
Y = (-radius(i,2):(radius(i,2)))' * ones(1,2*radius(i,1)+1);
Z =  X.^2/radius(i,1)^2 + Y.^2/radius(i,2)^2;
Znew = Z;
% gives each object its Object ID
Znew(Z <= 1) = CurrentModelMatrix.ObjectID(i);
Znew(Z > 1) = 0;

% but there might be already other objects in the cut out image...Keep them
excission = virtualimage(upperboundary(i):lowerboundary(i),leftboundary(i):rightboundary(i));
coordinates = excission == CurrentModelMatrix.ObjectID(i);
excission(coordinates) =0;
Znew = uint16(Znew);

 %idea: rotation only for nuclei in the middle of the image 
 %treatment of border case
if (leftboundary(i) > 1 && upperboundary(i) >1 && rightboundary(i)<columnnumber && lowerboundary(i)<rownumber)
   beforerotation=Znew;
   afterrotation=imtransform(beforerotation,tform);
   %some values are set by the rotation to other values than the objectID.
   %Correct this
   afterrotation(afterrotation > 0)=CurrentModelMatrix.ObjectID(i);
   dimensionofnewimage=size(afterrotation);

    if ((center(i,2)-ceil(dimensionofnewimage(1)/2)>=1) && (center(i,2)+floor(dimensionofnewimage(1)/2)-1)<=rownumber && center(i,1)-ceil(dimensionofnewimage(2)/2)>=1 && (center(i,1)+floor(dimensionofnewimage(2)/2)-1)<=columnnumber)
        virtualimage((center(i,2)-ceil(dimensionofnewimage(1)/2):center(i,2)+floor(dimensionofnewimage(1)/2)-1),(center(i,1)-ceil(dimensionofnewimage(2)/2):center(i,1)+floor(dimensionofnewimage(2)/2)-1))=virtualimage((center(i,2)-ceil(dimensionofnewimage(1)/2):center(i,2)+floor(dimensionofnewimage(1)/2)-1),(center(i,1)-ceil(dimensionofnewimage(2)/2):center(i,1)+floor(dimensionofnewimage(2)/2)-1))+ afterrotation;
    else
        virtualimage(upperboundary(i):lowerboundary(i), leftboundary(i):rightboundary(i)) = beforerotation;
    end
else
    if ~logical((leftboundary(i)-1)/(leftboundary(i)-1+eps) )
        Znew = Znew(:,(2*radius(i,1)+1)-(rightboundary(i)-leftboundary(i)):(2*radius(i,1)+1));
    end
    if ~logical((rightboundary(i)-columnnumber)/(leftboundary(i)-columnnumber+eps) )
        Znew = Znew(:, 1:(1+(rightboundary(i)-leftboundary(i))));
    end
    if ~logical((upperboundary(i)-1)/(upperboundary(i)-1+eps) )
       Znew =  Znew((2*radius(i,2)+1)-(lowerboundary(i)-upperboundary(i)):2*radius(i,2)+1, : );
    end
    if ~logical((lowerboundary(i)-rownumber)/(lowerboundary(i)-rownumber+eps) )
        Znew = Znew(1:(1+(lowerboundary(i)-upperboundary(i))),:);
    end
    virtualimage(upperboundary(i):lowerboundary(i), leftboundary(i):rightboundary(i))= excission + Znew;
end
end


