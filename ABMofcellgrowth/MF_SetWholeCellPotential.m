function[energymatrix]=MF_SetWholeCellPotential(CurrentModelMatrix,interactionmatrix, iCell,...
    occupied_area,nuclei_area, nuclei_outlines, energymatrix,normalized,Store)
%[MF] seems fine but I must come back on this one 06/13
%for each cell gives out the potential function of the distance away from
%the nucleus membrane, gives a kind of Heavyside function for nucleus and
%other cells.

center = round(CurrentModelMatrix.Nuclei_Location(iCell,:));
radius = CurrentModelMatrix.EnRadius;
columnnumber = CurrentModelMatrix.columnnumber;
rownumber = CurrentModelMatrix.rownumber;
%defines border values of image which is cut out 
leftboundary  = center(:,1)-radius;
leftboundary(leftboundary<1)=1;
newrad_1=center(:,1)-leftboundary;
rightboundary = center(:,1)+radius;
rightboundary(rightboundary>columnnumber)=columnnumber;
newrad_2=rightboundary-center(:,1);
upperboundary = center(:,2)-radius;
upperboundary(upperboundary<1)=1;
newrad_3=center(:,2)-upperboundary;
lowerboundary = center(:,2)+radius;
lowerboundary(lowerboundary>rownumber)=rownumber;
newrad_4=lowerboundary-center(:,2);

%field of view

%area of own cell (specific_area)
cutout2=nuclei_area(upperboundary:lowerboundary,leftboundary:rightboundary);

%surroundings
    cutoutNucl = double(cutout2);
    cutoutNucl(cutoutNucl == CurrentModelMatrix.ObjectID(iCell))=0;
    potentialheight_nucleus=1;
    cutoutNucl(cutoutNucl>0)=potentialheight_nucleus;
    energymatrix.pt_allsurroundingNucl=cutoutNucl;
%%%

specific_area=double(cutout2);
specific_area(specific_area~=CurrentModelMatrix.ObjectID(iCell))=0;
SE3 = strel('disk', 5);
specific_area=imdilate(specific_area,SE3);
potentialheight_nucleus=1;
specific_area(specific_area>0)=potentialheight_nucleus;
%figure;imshow(specific_area,[]);impixelinfo
energymatrix.pt_nuclei_area=specific_area;

%area of other cells
cutout3=occupied_area(upperboundary:lowerboundary,leftboundary:rightboundary);%double(
cutout3(cutout3==CurrentModelMatrix.ObjectID(iCell))=0;
%figure;imshow(cutout3,[]);impixelinfo

if Store == false
energymatrix.pt_other_cells=logical(cutout3 + cutoutNucl);
elseif Store == true
energymatrix.StoredCellpotential{iCell}=logical(cutout3 + cutoutNucl);
end


%calculate position of nodes for sum_of_potential calculations
%test:: nuclei_area=CurrentModelMatrix.nuclei_area;
bordernode_location=nuclei_area;
pos = find(bordernode_location > 0);
bordernode_location(pos)=200; %#ok<FNDSB>
nodesofcell=interactionmatrix.bordernodes{iCell,1};
nodesofcell(nodesofcell(:,5)>0,:)=[];

for i = 1:size(nodesofcell,1)
coordinate_condition3=nodesofcell(i,3);
coordinate_condition4=nodesofcell(i,4);
checkcoordinate_condition=~reshape((coordinate_condition3 > 1 & coordinate_condition3 < ...
    CurrentModelMatrix.columnnumber & coordinate_condition4 > 1 & coordinate_condition4 < CurrentModelMatrix.rownumber),[],1);
coordinatesfor1=reshape(CurrentModelMatrix.rownumber.*nodesofcell(i,3)+nodesofcell(i,4),[],1);
coordinatesfor1(checkcoordinate_condition)=[];
bordernode_location(coordinatesfor1)=10000+i; %here is this story of 400 ... that is baaaad let's change it
end

cutout5=bordernode_location(upperboundary:lowerboundary,leftboundary:rightboundary);
energymatrix.currentnodelocations=double(cutout5);

cutout=nuclei_outlines(upperboundary:lowerboundary,leftboundary:rightboundary);
specific_outlines=cutout;
specific_outlines(specific_outlines~=CurrentModelMatrix.ObjectID(iCell))=0;

    
SE8 = strel('disk', 8);
largeextendedimage=imdilate(specific_outlines,SE8);

[r,c]=size(largeextendedimage);
if newrad_1<radius
  largeextendedimage(:,1)=CurrentModelMatrix.ObjectID(iCell);
end

if newrad_2<radius
    largeextendedimage(:,c)=CurrentModelMatrix.ObjectID(iCell);
end

if newrad_3<radius
   largeextendedimage(1,:)=CurrentModelMatrix.ObjectID(iCell);
end

if newrad_4<radius
    largeextendedimage(r,:)=CurrentModelMatrix.ObjectID(iCell);
end
       
largeextendedimage=imfill(largeextendedimage,'holes');

largeextendedimage(:,1)=0;
largeextendedimage(:,c)=0;
largeextendedimage(1,:)=0;
largeextendedimage(r,:)=0;

outer_outline=edge(largeextendedimage,'canny');
outer_outline=double(outer_outline);
outer_outline(outer_outline>0)=CurrentModelMatrix.ObjectID(iCell);
  

distancefromnucleimembrane = double(bwdist(outer_outline));
%figure;surf(distancefromnucleimembrane)
%distancefromnucleimembrane(specific_area>0)=0;
%possiblity to flatten potentialfct if potentialfct value of cells are set
%inside those cells to one value...=> normalization
%distancefromnucleimembrane(cutout3>0)=0;
distancefromnucleimembrane_potential=distancefromnucleimembrane.^3;%also quite driven, let see later. > it is what limits the expension
% figure;surf(distancefromnucleimembrane_potential)



%Display purpose: make outlines appear...
%virtualimage=visualize_celloutline_model_cell(interactionmatrix,CurrentModelMatrix);
%cutout4=double(virtualimage(upperboundary(iCell):lowerboundary(iCell),leftboundary(iCell):rightboundary(iCell)));
%cutout4(cutout4~=300)=0;
%cutout4(cutout4==300)=0.3;

%figure;imshow(distancefromnucleimembrane,[])
 [temp1,temp2]=size(distancefromnucleimembrane_potential);
if normalized==1
    
%     if temp1>225 && temp2>225
%     normfactor_here=max(max(distancefromnucleimembrane_potential(25:225,25:225)));
%     else
    normfactor_here=max(max(distancefromnucleimembrane_potential));
%     end
    cell_distancepotential=(distancefromnucleimembrane_potential)/normfactor_here;
else
   cell_distancepotential=(distancefromnucleimembrane_potential); 
end

% BasicLandscape = CurrentModelMatrix.basicLandscape(1:temp1,1:temp2);
% BasicLandscape = BasicLandscape./3;
% 
% cell_distancepotential = BasicLandscape + cell_distancepotential;
 cell_distancepotential = cell_distancepotential./max(cell_distancepotential(:));
%Bigger = 10.*nan(250);Bigger(26:226,26:226) = cell_distancepotential;Bigger(Bigger>0.157 | isnan(Bigger)) = 0.157;
energymatrix.current_distancepotential=cell_distancepotential;


%cutout4(cell_distancepotential>0.4)=0;
%with_cell_outlines=cell_distancepotential+cutout4;
%figure;mesh(Bigger)
%figure; mesh(cell_distancepotential)
%figure; plot(cell_distancepotential(100,1:200))

%figure;surfc(with_cell_outlines)





