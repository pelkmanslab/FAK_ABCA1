function [list] = MF_FindNucleiTension_list(CurrentModelMatrix,Trshld,whichcell)%, virtualimage)
%[MF]Massive changes 06/13 to the bottom
%calculates the overlapping nucleis. Gives out a list with 
%[Nr. Nucleus 1, Nr. Nucleus 2, angle between Nucleus 1&2]
TrshldOverlap = 10;
if nargin == 3
%get all cells which are concerned by overlap 
listofcellstotal=[];
for iCells = whichcell
listofcells= MF_ListObjectsAroundObject(CurrentModelMatrix, iCells, 400);
listofcellstotal=[listofcellstotal;listofcells];
end
listofcellstotal=unique(listofcellstotal');
else
listofcellstotal = 1:CurrentModelMatrix.numberofcells;    
end
if nargin >= 2
 TrshldOverlap = Trshld;
end

overlap_distancematrix=zeros(CurrentModelMatrix.numberofcells, CurrentModelMatrix.numberofcells);
overlap_angle=zeros(CurrentModelMatrix.numberofcells, CurrentModelMatrix.numberofcells);

%Calculate the distance from the centroid of the nucleus to the nucleus
%membrane in the direction (angle) where the other cell is.

for nucleiWelookAt = listofcellstotal

listofnucleiinrange= MF_ListObjectsAroundObject(CurrentModelMatrix, nucleiWelookAt, 400);%400 as a threshold, deserves some improvement

angle=zeros(1,length(listofnucleiinrange));
arctanindex= (-CurrentModelMatrix.Nuclei_Location(nucleiWelookAt,1)+CurrentModelMatrix.Nuclei_Location(listofnucleiinrange,1))>0;
arctanindexCAL=listofnucleiinrange(arctanindex);%nuclei shifted towards inf on X...
angle(arctanindex)=atan((CurrentModelMatrix.Nuclei_Location(nucleiWelookAt,2)-CurrentModelMatrix.Nuclei_Location(arctanindexCAL,2))./(-CurrentModelMatrix.Nuclei_Location(nucleiWelookAt,1)+CurrentModelMatrix.Nuclei_Location(arctanindexCAL,1)));
arctanindexneg=find(-CurrentModelMatrix.Nuclei_Location(nucleiWelookAt,1)+CurrentModelMatrix.Nuclei_Location(listofnucleiinrange,1)<0);
arctanindexnegCAL=listofnucleiinrange(arctanindexneg);%nuclei shifted towards 0 on X...
arctanindexzero=find(-CurrentModelMatrix.Nuclei_Location(nucleiWelookAt,1)+CurrentModelMatrix.Nuclei_Location(listofnucleiinrange,1)==0);

%for those at the same level
for index_arctanindex=1:length(arctanindexzero)
if -CurrentModelMatrix.Nuclei_Location(nucleiWelookAt,2)+CurrentModelMatrix.Nuclei_Location(listofnucleiinrange(arctanindexzero(index_arctanindex)),2)>0
    angle(arctanindexzero(index_arctanindex))= pi * 3/2;
else
    angle(arctanindexzero(index_arctanindex)) = pi/2;
end
end

%angle calculation for the other cases
angle(arctanindexneg)=pi - atan((CurrentModelMatrix.Nuclei_Location(nucleiWelookAt,2)-CurrentModelMatrix.Nuclei_Location(arctanindexnegCAL,2))./(CurrentModelMatrix.Nuclei_Location(nucleiWelookAt,1)-CurrentModelMatrix.Nuclei_Location(arctanindexnegCAL,1)));
angle(angle<0)=2*pi+angle(angle<0);
%nice

%calculate the distance from centroid to nuclei membrane in each direction where there is a close cell   
factorl=zeros(length(angle),1);
centertonucleimembrane=zeros(1,length(angle));
for incre=1:length(angle)
    factort=atan(CurrentModelMatrix.radius(nucleiWelookAt,1)/CurrentModelMatrix.radius(nucleiWelookAt,2)*tan(angle(incre)-CurrentModelMatrix.angle(nucleiWelookAt,1)));
    factorl(incre,1)=sqrt((CurrentModelMatrix.radius(nucleiWelookAt,1)*cos(factort))^2 + (CurrentModelMatrix.radius(nucleiWelookAt,2)*sin(factort))^2);
    centertonucleimembrane(1,incre)=factorl(incre,1);
end   


if ~isempty(listofnucleiinrange)
for i=1:length(listofnucleiinrange)
    %there is options for optimization here.
overlap_distancematrix(nucleiWelookAt,listofnucleiinrange(i)) = centertonucleimembrane(i); % limit to define an overlap, huge matrix after a while
overlap_angle(nucleiWelookAt,listofnucleiinrange(i))=angle(i);%0 values represents cells that are too far to be considered 
end 

end 

end 
clear centertonucleimembrane angle



%Add now the distance of the two centroid-nuclei membrane together
%let's do something more efficient, and go just where it is worth it
posTMP = overlap_distancematrix ~= 0;
[IdxTMP(:,1) IdxTMP(:,2)] = find(posTMP);
overlap_distance = [zeros(sum(posTMP(:)),1) IdxTMP];
clear posTMP IdxTMP
%let's define the distance between the two centroids that must not be reached
for i = 1:size(overlap_distance,1)        
overlap_distance(i,1) = overlap_distancematrix(overlap_distance(i,2),overlap_distance(i,3))... %min dist to be in the nuclei
    +overlap_distancematrix(overlap_distance(i,3),overlap_distance(i,2)); %the same but for the considered close cell in that direction      
end

%comparison of the actual distance between centroids with the limit
%distance calculated just above

%init.
nucleidistance_matrix = overlap_distance;
nucleidistance_matrix(:,1) = 0;

%euclidian distance between centroids of interrest
for i = 1:size(nucleidistance_matrix,1)
nucleidistance_matrix(i,1)= ((CurrentModelMatrix.Nuclei_Location(nucleidistance_matrix(i,3),1)-CurrentModelMatrix.Nuclei_Location(nucleidistance_matrix(i,2),1)).^2 ...
    +(CurrentModelMatrix.Nuclei_Location(nucleidistance_matrix(i,3),2)-CurrentModelMatrix.Nuclei_Location(nucleidistance_matrix(i,2),2)).^2).^(1/2);
end

%Look at the difference
TMPdiffList =  nucleidistance_matrix(:,1) - (overlap_distance(:,1)+TrshldOverlap);%everything that is at less than 15 pixels than the treshold.
posTMP = TMPdiffList < 0;
list = nucleidistance_matrix(posTMP,2:3) ;
list = [list zeros(size(list,1),1)];
%than stock the angle that corresponds
for i = 1:size(list,1)
    list(i,3) = overlap_angle(list(i,1),list(i,2));
end
% clear TMPdiffList posTMP
% displTXT = size(list,1);
%msgTMP = sprintf('there is %d overlapping nuclei',displTXT);
% disp(msgTMP);




%plotting of nucleis which overlap

%if ~isempty(list)
% 10*sqrt((sin(2*diffangle-pi/4)).^2)
% test: list=tension_list;
%objectIDList(:,1)=input_matrix.ObjectID(list(:,1));
%objectIDList(:,2)=input_matrix.ObjectID(list(:,2));

%colortestimage1=virtualimage;
%for i=objectIDList'
%    colortestimage1(find(colortestimage1==i(1,1)))=400;
%    colortestimage1(find(colortestimage1==i(2,1)))=400;
%end

%colortestimage1(find(colortestimage1>0 & colortestimage1~=400))=200;
%[r,q]=size(colortestimage1);
%colortestimage2a=zeros(r,q);
%colortestimage2b=colortestimage2a;
%colortestimage2c=colortestimage2a;

%colortestimage2a(colortestimage1==400)=1;
%colortestimage2b(colortestimage1==400)=0.6;
%colortestimage2c(colortestimage1==400)=0.2;
 
%colortestimage2a(colortestimage1==200)=0.2;
%colortestimage2b(colortestimage1==200)=0.6;
%colortestimage2c(colortestimage1==200)=0.8;

%colorimage=cat(3,colortestimage2a,colortestimage2b,colortestimage2c);
%figure
%imshow(colorimage)   

%end
