function[CurrentModelMatrix,interactionmatrix,energymatrix]=MF_TOcheck(CurrentModelMatrix, interactionmatrix,energymatrix)
%[MF]shortly look at 06/13 inserted the FrozenTag
minimal_dist_to_nucleimembrane=3;

% iteration2=[];iteration1=[];

%what to do if more extremaList points...


 for iteration2=1:CurrentModelMatrix.numberofcells
 if ~CurrentModelMatrix.FreezeTag(iteration2)
 for iteration1=1:CurrentModelMatrix.numberofcells
 if ~CurrentModelMatrix.FreezeTag(iteration1)
         
if ~isempty(interactionmatrix.extremaList{iteration2,iteration1})       
 iCell=iteration2;



 
%iCell=13;
%iteration1=224;
newnodeNb=interactionmatrix.numberofbordernodes(iCell,1)+1;
newnodeNb2=interactionmatrix.numberofbordernodes(iCell,1)+2;
tempvalues=interactionmatrix.extremaList{iCell,iteration1};
%do not place the same extrema point again....
if find(interactionmatrix.bordernodes{iCell,1}(:,3)==round(tempvalues(1,1)))
    if find(interactionmatrix.bordernodes{iCell,1}(:,4)==round(tempvalues(1,2)))
        if find(interactionmatrix.bordernodes{iCell,1}(:,2)==tempvalues(1,3))
            continue
        end
    end
end
interactionmatrix.bordernodes{iCell,1}(newnodeNb,3)=round(tempvalues(1,1));
interactionmatrix.bordernodes{iCell,1}(newnodeNb,4)=round(tempvalues(1,2));
interactionmatrix.bordernodes{iCell,1}(newnodeNb,2)=tempvalues(1,3);
interactionmatrix.bordernodes{iCell,1}(newnodeNb,5)=2;
interactionmatrix.bordernodes{iCell,1}(newnodeNb,1)=sqrt((CurrentModelMatrix.Nuclei_Location(iCell,1)-round(tempvalues(1,1))).^2+(CurrentModelMatrix.Nuclei_Location(iCell,2)-round(tempvalues(1,2))).^2);

interactionmatrix.bordernodes{iCell,1}(newnodeNb2,3)=round(tempvalues(2,1));
interactionmatrix.bordernodes{iCell,1}(newnodeNb2,4)=round(tempvalues(2,2));
interactionmatrix.bordernodes{iCell,1}(newnodeNb2,2)=tempvalues(2,3);
interactionmatrix.bordernodes{iCell,1}(newnodeNb2,5)=2;
interactionmatrix.bordernodes{iCell,1}(newnodeNb2,1)=sqrt((CurrentModelMatrix.Nuclei_Location(iCell,1)-round(tempvalues(2,1))).^2+(CurrentModelMatrix.Nuclei_Location(iCell,2)-round(tempvalues(2,2))).^2);
%calculate here the minimal distance 
factorl=zeros(interactionmatrix.numberofnodes,1);
dist1=factorl;
dist2=factorl;
angles=interactionmatrix.bordernodes{iCell,1}(:,2);
for node=[newnodeNb,newnodeNb2]
    factort=atan(CurrentModelMatrix.radius(iCell,1)/CurrentModelMatrix.radius(iCell,2)*tan(angles(node)-CurrentModelMatrix.angle(iCell,1)));
    factorl(node,1)=sqrt((CurrentModelMatrix.radius(iCell,1).*cos(factort)).^2 + (CurrentModelMatrix.radius(iCell,2).*sin(factort)).^2);
    dist1(node,1)=factorl(node)*cos(angles(node));
    dist2(node,1)=factorl(node)*sin(angles(node));
    interactionmatrix.bordernodes{iCell}(node,6)=factorl(node,1)+minimal_dist_to_nucleimembrane;
    %interactionmatrix.centertonucleimembrane(node,1,i)=factorl(node,1);
end
%sort rows
interactionmatrix.bordernodes{iCell}=sortrows(interactionmatrix.bordernodes{iCell},2);
interactionmatrix.numberofbordernodes(iCell)=length(interactionmatrix.bordernodes{iCell,1}(:,2));
end
 end
 end
end
end
%iCell=10;
%Transfer of extremapoint complete, delete them
interactionmatrix.extremaList=cell(CurrentModelMatrix.numberofcells,CurrentModelMatrix.numberofcells);


for iCell=1:CurrentModelMatrix.numberofcells
 if ~CurrentModelMatrix.FreezeTag(iCell)
    %iCell=21
%find out which nodes have to be retracted to connecting line.
index2= interactionmatrix.bordernodes{iCell,1}(:,7)>0&interactionmatrix.bordernodes{iCell,1}(:,5)>1;
interactionmatrix.bordernodes{iCell,1}(index2,7)=0;
index1=find(interactionmatrix.bordernodes{iCell,1}(:,7)>0);

for i=1:length(index1)
    node=index1(i);
    
    %if there is an overlap, first check with the corresponding node, whether
    %it is an innernode, or if the node is already at minimal distance...
    %not_retractable=0;

    %if interactionmatrix.bordernodes{iCell,1}(node,5)==1
    %not_retractable=1;
    %end

    %if interactionmatrix.bordernodes{iCell,1}(node,1)<interactionmatrix.bordernodes{iCell,1}(node,6)+1
    %not_retractable=1;
    %end

    %if not_retractable==0
    lowerboundary=find(round(interactionmatrix.bordernodes{iCell,1}(1:node,5))>1, 1, 'last' );
    %lowerboundary=max(find(round(interactionmatrix.bordernodes{iCell,1}(1:node,5))==2));
  
    if isempty(lowerboundary)
        lowerboundary=find(round(interactionmatrix.bordernodes{iCell,1}(:,5))>1, 1, 'last' );
    end
    upperboundary=find(round(interactionmatrix.bordernodes{iCell,1}((node+1):end,5))>1, 1 )+node;
    if isempty(upperboundary)
        upperboundary=find(round(interactionmatrix.bordernodes{iCell,1}(:,5))>1, 1 );
    end
    
    p1x=interactionmatrix.bordernodes{iCell,1}(lowerboundary,3);
    p1y=interactionmatrix.bordernodes{iCell,1}(lowerboundary,4);
    p2x=interactionmatrix.bordernodes{iCell,1}(upperboundary,3);
    p2y=interactionmatrix.bordernodes{iCell,1}(upperboundary,4);
    c1=CurrentModelMatrix.Nuclei_Location(iCell,1);
    c2=CurrentModelMatrix.Nuclei_Location(iCell,2);
    n1=interactionmatrix.bordernodes{iCell,1}(node,3);
    n2=interactionmatrix.bordernodes{iCell,1}(node,4);
    
    kappa=((c1-p1x)*(p2y-p1y) - (c2-p1y)*(p2x-p1x))/((n2-c2)*(p2x-p1x)-(n1-c1)*(p2y-p1y));
    interactionmatrix.bordernodes{iCell,1}(node,1)=sqrt((CurrentModelMatrix.Nuclei_Location(iCell,1)-interactionmatrix.bordernodes{iCell,1}(node,3)).^2+(CurrentModelMatrix.Nuclei_Location(iCell,2)-interactionmatrix.bordernodes{iCell,1}(node,4)).^2);
     
    if interactionmatrix.bordernodes{iCell,1}(node,1)*kappa>interactionmatrix.bordernodes{iCell,1}(node,6)
    
        factorl=interactionmatrix.bordernodes{iCell,1}(node,1)*kappa;
        interactionmatrix.bordernodes{iCell,1}(node,1)=factorl;
        interactionmatrix.bordernodes{iCell,1}(node,3)=round(factorl*cos(interactionmatrix.bordernodes{iCell,1}(node,2))+CurrentModelMatrix.Nuclei_Location(iCell,1));
        interactionmatrix.bordernodes{iCell,1}(node,4)=round(-factorl*sin(interactionmatrix.bordernodes{iCell,1}(node,2))+CurrentModelMatrix.Nuclei_Location(iCell,2));          
    else
        %reduce node to minimal length before you change other nodes...
        if interactionmatrix.bordernodes{iCell,1}(node,5)==1 %do not reduce length of innernodes
            factorl=interactionmatrix.bordernodes{iCell,1}(node,1);
        else
        factorl=interactionmatrix.bordernodes{iCell,1}(node,6);
        end
        interactionmatrix.bordernodes{iCell,1}(node,1)=factorl;
        interactionmatrix.bordernodes{iCell,1}(node,3)=round(factorl*cos(interactionmatrix.bordernodes{iCell,1}(node,2))+CurrentModelMatrix.Nuclei_Location(iCell,1));
        interactionmatrix.bordernodes{iCell,1}(node,4)=round(-factorl*sin(interactionmatrix.bordernodes{iCell,1}(node,2))+CurrentModelMatrix.Nuclei_Location(iCell,2));    
        %change now nodes of other cell
        othernuclei=interactionmatrix.bordernodes{iCell,1}(node,7);
        newnode=length(interactionmatrix.bordernodes{othernuclei,1}(:,2))+1;
        interactionmatrix.numberofbordernodes(othernuclei,1)=newnode;
        factorl=interactionmatrix.bordernodes{iCell,1}(node,1)*1.05;
        interactionmatrix.bordernodes{othernuclei,1}(newnode,3)=round(factorl*cos(interactionmatrix.bordernodes{iCell,1}(node,2))+CurrentModelMatrix.Nuclei_Location(iCell,1));
        interactionmatrix.bordernodes{othernuclei,1}(newnode,4)=round(-factorl*sin(interactionmatrix.bordernodes{iCell,1}(node,2))+CurrentModelMatrix.Nuclei_Location(iCell,2));
        interactionmatrix.bordernodes{othernuclei,1}(newnode,1)=sqrt((CurrentModelMatrix.Nuclei_Location(iCell,1)-interactionmatrix.bordernodes{othernuclei,1}(newnode,3)).^2+(CurrentModelMatrix.Nuclei_Location(iCell,2)-interactionmatrix.bordernodes{othernuclei,1}(newnode,4)).^2);
        coordinate2=interactionmatrix.bordernodes{othernuclei,1}(newnode,3);
        coordinate1=interactionmatrix.bordernodes{othernuclei,1}(newnode,4);
        %update now interactionmatrix for all nodes:
        %calculate angle
        if -CurrentModelMatrix.Nuclei_Location(othernuclei,1)+coordinate2>0
        angles=atan((CurrentModelMatrix.Nuclei_Location(othernuclei,2)-coordinate1)./(-CurrentModelMatrix.Nuclei_Location(othernuclei,1)+coordinate2));
        elseif -CurrentModelMatrix.Nuclei_Location(othernuclei,1)+coordinate2<0
        angles=pi- atan((CurrentModelMatrix.Nuclei_Location(othernuclei,2)-coordinate1)./(CurrentModelMatrix.Nuclei_Location(othernuclei,1)-coordinate2));
        else
            if -CurrentModelMatrix.Nuclei_Location(othernuclei,2)+coordinate1>0
            angles= pi * 3/2;
            else
            angles= pi/2;
            end
        end
        angles(angles<0)=2*pi+angles(angles<0);
        interactionmatrix.bordernodes{othernuclei,1}(newnode,2)=angles;
        
        %calculate here the minimal distance 
        factort=atan(CurrentModelMatrix.radius(othernuclei,1)/CurrentModelMatrix.radius(othernuclei,2)*tan(angles-CurrentModelMatrix.angle(othernuclei,1)));
        factorl=sqrt((CurrentModelMatrix.radius(othernuclei,1).*cos(factort)).^2 + (CurrentModelMatrix.radius(othernuclei,2).*sin(factort)).^2);
%         dist1=factorl*cos(angles);
%         dist2=factorl*sin(angles);
        interactionmatrix.bordernodes{othernuclei}(newnode,6)=factorl+minimal_dist_to_nucleimembrane;
        
        interactionmatrix.bordernodes{othernuclei}(newnode,5)=3;
        
       
            %replace 2 by 3 in column since edges are now important...
            lowerboundary=find(interactionmatrix.bordernodes{othernuclei,1}(:,2)<angles & round(interactionmatrix.bordernodes{othernuclei,1}(:,5))==2, 1, 'last' );
            if isempty(lowerboundary)
            lowerboundary=find(round(interactionmatrix.bordernodes{othernuclei,1}(:,5))==2, 1, 'last' );
            end
            upperboundary=find(interactionmatrix.bordernodes{othernuclei,1}(:,2)>angles &round(interactionmatrix.bordernodes{othernuclei,1}(:,5))==2, 1 );
            if isempty(upperboundary)
            upperboundary=find(round(interactionmatrix.bordernodes{othernuclei,1}(:,5))==2, 1 );
            end
            interactionmatrix.bordernodes{othernuclei,1}(lowerboundary,5)=2.1;
            interactionmatrix.bordernodes{othernuclei,1}(upperboundary,5)=2.1;
           
        if interactionmatrix.bordernodes{othernuclei}(newnode,1)<interactionmatrix.bordernodes{othernuclei,1}(newnode,6)
            %here move cells...
              move_vect=CurrentModelMatrix.Nuclei_Location(othernuclei,:)-CurrentModelMatrix.Nuclei_Location(iCell,:);
              direction_vect=-round(2*move_vect./(mean(abs(move_vect))));
    
              CurrentModelMatrix.Nuclei_Location(iCell,:)=CurrentModelMatrix.Nuclei_Location(iCell,:)+direction_vect;
              interactionmatrix.bordernodes{iCell}(:,3)=interactionmatrix.bordernodes{iCell}(:,3) + direction_vect(1);
              interactionmatrix.bordernodes{iCell}(:,4)=interactionmatrix.bordernodes{iCell}(:,4) + direction_vect(2);
              %check whether the nodes of the moved cell have to be placed
              %differently....(not to move cells into other cells...)
              %[overlaplist_temp, interactionmatrix]=Alternative_overlap_surrnodesSTD(CurrentModelMatrix, interactionmatrix);
              %if ~isempty(overlaplist_temp)
              %overlaplist_temp=overlaplist_temp(overlaplist_temp(:,2)==othernuclei,:);
              %[interactionmatrix,CurrentModelMatrix,energymatrix]=move_overlapping_nodes_againSTD(energymatrix,CurrentModelMatrix,interactionmatrix, overlaplist_temp);
              %end
        end
            
         %sort rows
        interactionmatrix.bordernodes{othernuclei}=sortrows(interactionmatrix.bordernodes{othernuclei},2);
        
    end %over a specific node
        
end %over all nodes of a cell
end
end  %over all cells



%delete all nodes introduced to sustain 
for i=1:CurrentModelMatrix.numberofcells
    if ~CurrentModelMatrix.FreezeTag(i)
  interactionmatrix.bordernodes{i,1}((interactionmatrix.bordernodes{i,1}(:,5)==2.1),5)=3;
  %interactionmatrix.bordernodes{i,1}((interactionmatrix.bordernodes{i,1}(:,5)==2),:)=[];
  %interactionmatrix.numberofbordernodes(i,1)=length(interactionmatrix.bordernodes{i,1}(:,2))
    end
end




