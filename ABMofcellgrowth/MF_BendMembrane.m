function [output1,interactionmatrix] = MF_BendMembrane(interactionmatrix,CurrentModelMatrix, iCell)
%inner nodes: place inner nodes so that no border lines cross nuclei!
%[MF] massive changes 06/13, what stephan call inner nodes refers to new
%nodes placed regularly around the nuclei of interrest, this is done when
%we need to recreate or just create new nodes that do generate lines
%crossing the nuclei.
output1=true;
if ~CurrentModelMatrix.FreezeTag(iCell)
output1=false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate the position of the innernodes...

%define a set of inner nodes.
cutStep = 18;%18 point to divide the circle, this value is crucial
incrementstep=2*pi/cutStep;
angleinnernode = 0:incrementstep:(2*pi - incrementstep); %hardcoded!
howmanyinnernodes=length(angleinnernode);

%calculate here how large is the distance to nuclei membrane from centroid  
factorl=zeros(howmanyinnernodes,1);
innernodeXpos=zeros(howmanyinnernodes,1);
innernodeYpos=zeros(howmanyinnernodes,1);

%defines positions of the innernodes, reminder for those who needs it: if
%A=tan(B) then B=atan(A). You are welcome ;)
for innernode=1:howmanyinnernodes
    factort=atan(CurrentModelMatrix.radius(iCell,1)/CurrentModelMatrix.radius(iCell,2) ...
        *tan(angleinnernode(innernode)-CurrentModelMatrix.angle(iCell,1)));
    factorl(innernode,1)=sqrt((CurrentModelMatrix.radius(iCell,1)*cos(factort))^2 + (CurrentModelMatrix.radius(iCell,2)*sin(factort))^2);
    factorl(innernode,1)=factorl(innernode,1)+5;%hardcoded too !!!
    dist1=factorl(innernode)*cos(angleinnernode(innernode));
    dist2=factorl(innernode)*sin(angleinnernode(innernode));
innernodeXpos(innernode,1)=round(CurrentModelMatrix.Nuclei_Location(iCell,1)+dist1);
innernodeYpos(innernode,1)=round(CurrentModelMatrix.Nuclei_Location(iCell,2)-dist2);
end   
%we now have all coordinates and angles and norms of the new nodes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if we look at a border cell, no need to waste time.
if ~isempty(find(interactionmatrix.bordernodes{iCell}<0, 1))
    output1=true;
%     output2=interactionmatrix ;
    return
end

%or if some inner node is negative:
if ~isempty(find(innernodeYpos<=0, 1))%.....here it was. -_-
    output1=true;
%     output2=interactionmatrix ;
    return
end

if ~isempty(find(innernodeXpos<=0, 1))%.....here it was. -_-
    output1=true;
%     output2=interactionmatrix ;
    return
end


%all angles of the current nodes
angle_nodesofcell=interactionmatrix.bordernodes{iCell}(:,2);

%is there an existing node already placed on an inner node?
AlignmentENonIN = false(howmanyinnernodes,1);
for innernodecheck=1:howmanyinnernodes
    temp1=angleinnernode(innernodecheck);%the angle of the IN of interrest
if ~isempty(find(temp1==angle_nodesofcell, 1))%is there a match
AlignmentENonIN(innernodecheck)=true;%yes so put a 1
end
end

%use a cell to control that specific inner nodes are only once checked!
%come back on this
% placingnodeList=cell(howmanyinnernodes,1);
% for i=1:howmanyinnernodes
%     placingnodeList{i}=i;
% end

placingnodeList=false(howmanyinnernodes,1);%false if not yet treated, true if treated


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while ~output1
%output control
  
lambda = ones(howmanyinnernodes,1);
%kabba= zeros(howmanyinnernodes,1);

for innernodecheck = 1:howmanyinnernodes

%do not calculate the values again for already placed nodes: trick use
%lambda(innernodecheck,1)=2 since those nodes are then not selected later in
%the model
if placingnodeList(innernodecheck)
    lambda(innernodecheck,1)=2;
    continue
end
   
%special case: innernode and placed node from interactionmatrix on same line    
if AlignmentENonIN(innernodecheck)
    lambda(innernodecheck,1)= 1;
    %kabba(innernodecheck,1)=0;
    continue
end
%at this point we are sure that nodes we will check are not aligned with
%existing node and are not checked a second time


%Calculate angles which are smaller or larger than the innernode to
%calculate the intersection point (then determine later
%if this node is in the cell)
if angleinnernode(innernodecheck)>interactionmatrix.bordernodes{iCell}(1,2)%thing to not forget, angles are stored from the smallest to the biggest
  firstnodetoconsider=sum(angle_nodesofcell<=angleinnernode(innernodecheck));%from all the nodes that are smaler, take the biggest
    if(firstnodetoconsider+1>length(angle_nodesofcell))%if the first node to consider is the last one (ie, there is no bigger angle)
    secondnodetoconsider=1;%take the first position
    else
    secondnodetoconsider= firstnodetoconsider+1;%you then have the position of the closest smaler and bigger
    end
else%if there is only angles bigger than angleinnernode(innernodecheck)
    firstnodetoconsider=interactionmatrix.numberofbordernodes(iCell);
    secondnodetoconsider=1;
end
%nodestoconsider=[firstnodetoconsider,secondnodetoconsider];

%get all necessary positions
p1x=CurrentModelMatrix.Nuclei_Location(iCell,1);
p1y=CurrentModelMatrix.Nuclei_Location(iCell,2);
p2x=innernodeXpos(innernodecheck);
p2y=innernodeYpos(innernodecheck);
p3x=interactionmatrix.bordernodes{iCell}(firstnodetoconsider,3);
p3y=interactionmatrix.bordernodes{iCell}(firstnodetoconsider,4);
p4x=interactionmatrix.bordernodes{iCell}(secondnodetoconsider,3);
p4y=interactionmatrix.bordernodes{iCell}(secondnodetoconsider,4);
%


if ((p2x-p1x)*(p4y-p3y)-(p2y-p1y)*(p4x-p3x))~=0 %if they are parallel/no division by zero!
 lambda(innernodecheck,1)= ((p1y - p3y)*(p4x-p3x)-(p1x-p3x)*(p4y-p3y))/((p2x-p1x)*(p4y-p3y)-(p2y-p1y)*(p4x-p3x));
else
 lambda(innernodecheck,1)=0;%of course lambda is not zero, but effectively, a node has to be placed...
end
%kabba(innernodecheck,1)=((p1x-p3x)+lambda(innernodecheck,1)*(p2x-p1x))/(p4x-p3x);
end

lambdatest=find(lambda<1);
innernodesequence=lambdatest;
%kappatest=find(kabba>1.3 | kabba<-0.3);
%innernodesequence=union(lambdatest,kappatest);
%if there are no "problem cases" terminate the module
if isempty(innernodesequence)
output1=true;
% output2=interactionmatrix ;
return
else
output1=false;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if there are lines which seem to cut the cells, here set appropriately new
%nodes.
%analysis of sequence. Find middle ones:1. make a difference list

%make a list of the difference between successive nodes keep in mind that the list is a circle
%exmpl: list : 15 16 18 24 difflist: 1 2 6 27 (there is 36 nodes defined at the begging, 0 does not count)

diffinnernodes=zeros(length(innernodesequence),1);
for i=1:length(innernodesequence)
    if i~=length(innernodesequence)
diffinnernodes(i)=innernodesequence(i+1)-innernodesequence(i);
    else
diffinnernodes(length(innernodesequence))=innernodesequence(1)-(innernodesequence(i)-cutStep);
    end
end

%Here the continuous regions of nodes are found. 
while ~isempty(innernodesequence)
    %first check whether 1 and 36 are both in the list, if so, you have to
    %be careful, since then "circular chains can occur and it is not
    %straigth forward
    
    if (innernodesequence(1)==1 &&innernodesequence(length(innernodesequence))==cutStep)
       ixx=find(diffinnernodes>=1,1,'last');%calculate the last element
       addnb=length(innernodesequence)-ixx;%how many elements belong to "circular chain from end"
       
       ix=find(diffinnernodes>=1,1,'first'); %find the first difference
       middleelement=(ix+addnb)/2; %how many elments are there in total / divided by two
        %check now whether the "middle" is at the end or at the beginning of the linear chain 
          if addnb>middleelement %are there more elements at the end or the beginning
              chosenodenb=ixx+middleelement;
          else
              chosenodenb=ix-middleelement;
          end
    else
        %no circular chain. Just determine from the beginning how many
        %nodes belong together.
        ix=find(diffinnernodes>=1,1,'first');
        middleelement=(ix)/2;
        chosenodenb=ix-middleelement;
    end
    
    %it could happen that there are equal nb of elements on both sites of
    %the angle
    chosenodenb=round(chosenodenb); %round to get position 
 %treat special case that chosenodbenb=0 meaning there are as many cells at the end as at the beginning.
     if chosenodenb==0
        chosenodenb=1;
      end
      
      chosenodenb=innernodesequence(chosenodenb);
      
      %now delete all node positions which have been considered above
      if (innernodesequence(1)==1 && innernodesequence(length(innernodesequence))==cutStep)
      diffinnernodes((ixx+1):length(innernodesequence))=[];
      innernodesequence((ixx+1):length(innernodesequence))=[];
      end
      
      diffinnernodes(1:ix)=[];
      innernodesequence(1:ix)=[];
   %save now all information of the newly added inner node. 
   placingnodeList(chosenodenb)=true;
   newnb=length(interactionmatrix.bordernodes{iCell,1}(:,2))+1;
   interactionmatrix.numberofbordernodes(iCell,1)=  newnb;
   extendfactor = factorl(chosenodenb,1) + 3;
   
   interactionmatrix.bordernodes{iCell}(newnb,1)= extendfactor; 
   dist1=extendfactor*cos(angleinnernode(chosenodenb));
   dist2=extendfactor*sin(angleinnernode(chosenodenb));
   extendednodeXpos=round(CurrentModelMatrix.Nuclei_Location(iCell,1)+dist1);
   extendednodeYpos=round(CurrentModelMatrix.Nuclei_Location(iCell,2)-dist2);
   
   interactionmatrix.bordernodes{iCell}(newnb,2)=angleinnernode(chosenodenb);
   interactionmatrix.bordernodes{iCell}(newnb,3)=extendednodeXpos;
   interactionmatrix.bordernodes{iCell}(newnb,4)=extendednodeYpos;
   interactionmatrix.bordernodes{iCell}(newnb,5)=1;
   interactionmatrix.bordernodes{iCell}(newnb,6)=extendfactor - 2;
   interactionmatrix.bordernodes{iCell}=sortrows(interactionmatrix.bordernodes{iCell},2);
%update now angle_nodesofcell for calculation (several innernodes might be
%placed per cell in this loop 
   angle_nodesofcell=interactionmatrix.bordernodes{iCell}(:,2);
   
end

 end
 %visualize here the outlines of the cell
 %virtualimage=visualize_celloutline_model_cell(interactionmatrix,CurrentModelMatrix);
end   
   

