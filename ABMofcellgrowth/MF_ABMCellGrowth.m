%Version June 2013 MF.
%matlabpool open 5
StartNumOfCells = 40;
NoN = 15;
interval = 0.25;
radius = 150;%can be considered as the maximal cell diameter.
maxCellsize = 2000;%nucleus area
minCellsize = 950;
%distribOfEnergy = [0.8,0.1];
totalenergy = [0.18,0.001]*NoN;%
nodeGrowthlimit = 60; %maximal distance a node can do in 1 hour. this value actually gives no constrain
foldername=('C:\Users\MAT\Desktop\ABM_newRound_Nature');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[virtualimage,CurrentModelMatrix,interactionmatrix,energymatrix] = MF_InitABmodel(StartNumOfCells, interval,totalenergy,NoN,radius);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[imagefolder,CurrentModelMatrixfolder,interactionmatrixfolder,energymatrixfolder,imagefolder2]=MF_SetFolders(foldername);

%%%%%%%%%%%%%%%%%%%
% % 
    CurrentModelMatrix= load('C:\Users\MAT\Desktop\ABM_newRound_Nature\CurrentModelMatrices\CurrentModelMatrix_95.png.mat');
    interactionmatrix= load('C:\Users\MAT\Desktop\ABM_newRound_Nature\interactionMatrices\interactionmatrix95.png.mat');
    energymatrix= load('C:\Users\MAT\Desktop\ABM_newRound_Nature\energymatrices\energymatrix95.png.mat');
    virtualimage = visualize_celloutline_model_cell(interactionmatrix,CurrentModelMatrix);
    foldername=('C:\Users\MAT\Desktop\ABM_newRound_Nature');
    [imagefolder,CurrentModelMatrixfolder,interactionmatrixfolder,energymatrixfolder,imagefolder2]=MF_SetFolders(foldername);

 for iteration_nb=95:600
     
     %the energy increases before being spent during the move, we do not
     %increase it after a certain value
     %if iteration_nb > 1
     pos = energymatrix.totalenergy(:,1) <= (0.5*NoN);
     energymatrix.totalenergy(pos) = energymatrix.totalenergy(pos,1) + (0.01*NoN*interval);
     %end
%if iteration_nb~=3    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%simulate cell growth (1 = interphase, 2 = mitotic cells)%%%%%%%%%

Bck1 = energymatrix;
Bck2 = CurrentModelMatrix;
Bck3 = interactionmatrix;

ReRun = true;
while ReRun == true
    
   energymatrix =  Bck1;
   CurrentModelMatrix =  Bck2;
   interactionmatrix =  Bck3;
   
    

disp('cells are about to grow')
tic
[CurrentModelMatrix, interactionmatrix]=MF_CellGrowthControl(energymatrix, CurrentModelMatrix, interactionmatrix, maxCellsize,minCellsize);
toc

disp('fitting node movement')
[interactionmatrix,CurrentModelMatrix,energymatrix,ReRun] = MF_SyncMovementOfNodes(energymatrix,CurrentModelMatrix,interactionmatrix,nodeGrowthlimit,inf);

end

% figure(2);imshow(IMGwithOverlap)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%simulate cell divsions
disp('cells are dividing')
mitoticcells=find(CurrentModelMatrix.mitotictag~=0);
interphasecells=find(CurrentModelMatrix.mitotictag==0);  
[readytodivide,CurrentModelMatrix] = MF_WhichCellMustDivide(CurrentModelMatrix);%, virtualimage); %this code define if a cell must divide
isx=intersect(interphasecells,readytodivide);



Bck1 = energymatrix;
Bck2 = CurrentModelMatrix;
Bck3 = interactionmatrix;
Bck4 = virtualimage;
Bck5 = isx;

ReRun = true;
while ReRun == true
    
   energymatrix =  Bck1;
   CurrentModelMatrix =  Bck2;
   interactionmatrix =  Bck3;
   virtualimage =  Bck4;
   isx =  Bck5;
   
tic
[CurrentModelMatrix,interactionmatrix,energymatrix]=MF_CellDivision(CurrentModelMatrix,interactionmatrix,energymatrix,virtualimage,isx,totalenergy,minCellsize);
toc

disp('trying to fit node movements coming from cell division')
[interactionmatrix,CurrentModelMatrix,energymatrix,ReRun] = MF_SyncMovementOfNodes(energymatrix,CurrentModelMatrix,interactionmatrix,nodeGrowthlimit,5);
end
% figure(3);imshow(IMGwithOverlap)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%move cells, here a constrain must be inserted about the movement allowed
%at each TP

Bck1 = energymatrix;
Bck2 = CurrentModelMatrix;
Bck3 = interactionmatrix;

ReRun = true;
while ReRun == true
    
   energymatrix =  Bck1;
   CurrentModelMatrix =  Bck2;
   interactionmatrix =  Bck3;
   
  
disp('the cells are moving, sensing the environment')
tic
[interactionmatrix,CurrentModelMatrix,energymatrix]=MF_SimCellMovement(energymatrix,CurrentModelMatrix,interactionmatrix,nodeGrowthlimit);
toc

disp('trying to fit node movements coming from spreading')
[interactionmatrix,CurrentModelMatrix,energymatrix,ReRun] = MF_SyncMovementOfNodes(energymatrix,CurrentModelMatrix,interactionmatrix,nodeGrowthlimit,50);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('all nodes correctly placed')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Update diverse parameters and status of cells
%calculate new LCD values
CurrentModelMatrix=MF_LCDinModel(CurrentModelMatrix);
    
%update status info of cells and mitotic tag info
CurrentModelMatrix.status = CurrentModelMatrix.status + CurrentModelMatrix.growthtimestepincrement;
CurrentModelMatrix.mitotictag(CurrentModelMatrix.mitotictag~=0)=CurrentModelMatrix.mitotictag(CurrentModelMatrix.mitotictag~=0)-1; 

for i=1:CurrentModelMatrix.numberofcells
interactionmatrix.bordernodes{i,1}(:,7)=zeros(interactionmatrix.numberofbordernodes(i),1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save and plot image
CurrentModelMatrixname=sprintf('CurrentModelMatrix_%d.png',iteration_nb);
interactionmatrixname=sprintf('interactionmatrix%d.png',iteration_nb);
energymatrixname=sprintf('energymatrix%d.png',iteration_nb);
save([CurrentModelMatrixfolder '\\' CurrentModelMatrixname '.mat'], '-struct','CurrentModelMatrix')
save([interactionmatrixfolder '\\' interactionmatrixname '.mat'],'-struct','interactionmatrix')
save([energymatrixfolder '\\' energymatrixname '.mat'],'-struct','energymatrix')

virtualimage=MF_ModelBckbone_Img(interactionmatrix,CurrentModelMatrix);
    
image2=zeros(CurrentModelMatrix.rownumber,CurrentModelMatrix.columnnumber);
image2(virtualimage>0)=300;
imagename=sprintf('actualimage%d.jpg',iteration_nb);
filename=[imagefolder '\\' imagename];
filename2=[imagefolder2 '\\' imagename];

imwrite(uint8(virtualimage),filename,'jpeg')
imwrite(uint8(image2),filename2,'jpeg')

 end

 
 
 
 
 
 
 
 
 
 
 

% % % [overlaplist, interactionmatrix,IMGwithOverlap]=new_idea_plot_overlap_more_detail(CurrentModelMatrix, interactionmatrix);% go again on it, check out angle definitions
% % % % imshow(IMGwithOverlap)
% % % %end
% % % while ~isempty(overlaplist)
% % % % if countReleaseAttempts < 200 ;
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%THis is all about node correct placement%%%%%%%%
% % % %First try to replace nodes "for easy cases". Use Standard procedure for
% % % %node placement (no nucleus movement)... 
% % % disp('first node displacement')
% % % if ~isempty(overlaplist)
% % %     [interactionmatrix,CurrentModelMatrix,energymatrix, worked]=...
% % %         move_overlapping_nodes_again_FAST_STD(energymatrix,CurrentModelMatrix,interactionmatrix, overlaplist,nodeGrowthlimit);
% % % end
% % % disp('second node displacement with nuclei shift')
% % % % %secondly deal with cases where nodes such as innernodes are involved 
% % % [overlaplist, interactionmatrix]=new_idea_plot_overlap_more_detail(CurrentModelMatrix, interactionmatrix);%, energymatrix);
% % % if ~isempty(overlaplist)
% % %     [CurrentModelMatrix, interactionmatrix,energymatrix]=...
% % %         correct_overlap_introduction_lineSTD(CurrentModelMatrix, interactionmatrix,energymatrix);
% % % end
% % % [overlaplist, interactionmatrix]=new_idea_plot_overlap_more_detail(CurrentModelMatrix, interactionmatrix);%, energymatrix);
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % % 
% % % count=0;
% % % while ~isempty(overlaplist) && count < 10
% % %     disp('last attempt before replacing nodes')
% % %     [interactionmatrix,CurrentModelMatrix,energymatrix,worked]=move_overlapping_nodes_again_FAST_STD(energymatrix,CurrentModelMatrix,interactionmatrix, overlaplist);
% % % if worked==0
% % %     count=count+1;
% % % end
% % % [overlaplist, interactionmatrix]=new_idea_plot_overlap_more_detail(CurrentModelMatrix, interactionmatrix);%, energymatrix);
% % % if ~isempty(overlaplist)
% % %     tic
% % %     [CurrentModelMatrix, interactionmatrix,energymatrix]=correct_overlap_introduction_lineSTD(CurrentModelMatrix, interactionmatrix,energymatrix);
% % %     toc
% % % end
% % % 
% % % [overlaplist, interactionmatrix]=new_idea_plot_overlap_more_detail(CurrentModelMatrix, interactionmatrix);%, energymatrix);
% % % count=count+1;
% % % disp('loop completed to reduce overlap')
% % % end
% % %  
% % % %it did not work so far? Replace all nodes of all cells involved and try again
% % % if ~isempty(overlaplist)
% % %     disp('fix hard cases with new node placing')
% % %     tic
% % % [interactionmatrix,CurrentModelMatrix,energymatrix]=model_move_nodes_selectionSTD(energymatrix,CurrentModelMatrix,interactionmatrix,nodeGrowthlimit);
% % % toc
% % % [overlaplist, interactionmatrix,IMGwithOverlap]=new_idea_plot_overlap_more_detail(CurrentModelMatrix, interactionmatrix);%, energymatrix);
% % % end
% % % 
% % %  
% % % % else
% % % % %maximalcontraction
% % % % MaxReleaseList = unique(overlaplist(:));
% % % %     for iCell = MaxReleaseList';
% % % %     [interactionmatrix] = init_nodes(interactionmatrix,CurrentModelMatrix, iCell,1);
% % % %     end
% % % %     for regrow = 1:50
% % % %     [interactionmatrix,CurrentModelMatrix,energymatrix]=model_move_nodesSTD(energymatrix,CurrentModelMatrix,interactionmatrix);
% % % %     end
% % % %     [overlaplist, interactionmatrix]=new_idea_plot_overlap_more_detail(CurrentModelMatrix, interactionmatrix);
% % % %     countReleaseAttempts = 1;
% % % % end
% % % % [interactionmatrix,CurrentModelMatrix]=NonIterativeCrashFix(IMGwithOverlap,CurrentModelMatrix,interactionmatrix);
% % % % countReleaseAttempts = countReleaseAttempts +1;
% % % % [overlaplist, interactionmatrix,IMGwithOverlap]=new_idea_plot_overlap_more_detail(CurrentModelMatrix, interactionmatrix);
% % % % figure(countReleaseAttempts);imshow(IMGwithOverlap)
% % % end





% % % %release  the population tensions, in order to have
% % % %growing islets.
% % % 
% % % disp('Now relax islets tensions')
% % % Bck1 = CurrentModelMatrix;
% % % Bck2 = interactionmatrix;
% % % ReRun = true;
% % % while ReRun == true
% % %     
% % %    CurrentModelMatrix =  Bck1;
% % %    interactionmatrix  =  Bck2;
% % %    ReRun = false;
% % % overlapping_list= get_tensions_in_islets(CurrentModelMatrix);%, virtualimage2
% % % %proceed with tension release as long as there are still overlapping nuclei
% % % while ~isempty(overlapping_list)
% % % %release
% % % CurrentModelMatrix = release_Islet_tensions(CurrentModelMatrix,interactionmatrix, overlapping_list);
% % % %determine which nuclei overlap
% % % overlapping_list= get_tensions_in_islets(CurrentModelMatrix);%, virtualimage2)
% % % end
% % % 
% % % disp('Fixing resulting overlaps')
% % % [overlaplist, interactionmatrix]=new_idea_plot_overlap_more_detail(CurrentModelMatrix, interactionmatrix);
% % %  %figure(2);imshow(IMGwithOverlap)
% % % %second release
% % % count = 0;
% % % while ~isempty(overlaplist)
% % %     count = count+1;test = (1:1000)';
% % % %      if any((count./5) == test)
% % % %          param = [0,1];
% % % % %      elseif any((count./10)==test)
% % % %      else
% % %          param = [0,0];
% % % %      end
% % %      
% % %      
% % %     if  count == 20
% % %         ReRun = true;
% % %         disp('redo Islet tension release')
% % %          break
% % %     end
% % %     
% % %     tic
% % %     [interactionmatrix,CurrentModelMatrix,energymatrix, worked]=...
% % %         move_overlapping_nodes_again_FAST_STD(param,energymatrix,CurrentModelMatrix,interactionmatrix, overlaplist,nodeGrowthlimit);
% % %     [overlaplist, interactionmatrix]=new_idea_plot_overlap_more_detail(CurrentModelMatrix, interactionmatrix);
% % %     %figure(1);imshow(IMGwithOverlap)
% % %     toc 
% % %     
% % %     
% % %     if ~isempty(overlaplist)%any((count./10) == test)
% % %     [CurrentModelMatrix, interactionmatrix,energymatrix]=correct_overlap_introduction_lineSTD(CurrentModelMatrix, interactionmatrix,energymatrix);  
% % %     [overlaplist, interactionmatrix,IMGwithOverlap]=new_idea_plot_overlap_more_detail(CurrentModelMatrix, interactionmatrix);
% % %     %figure(1);imshow(IMGwithOverlap)
% % %     end
% % %     
% % %     
% % %     if ~isempty(overlaplist)
% % %     [interactionmatrix,CurrentModelMatrix,energymatrix]=model_move_nodes_selectionSTD(energymatrix,CurrentModelMatrix,interactionmatrix,nodeGrowthlimit);
% % %     [overlaplist, interactionmatrix]=new_idea_plot_overlap_more_detail(CurrentModelMatrix, interactionmatrix);
% % %     end
% % % end
% % % 
% % % end
% % % 
















 
