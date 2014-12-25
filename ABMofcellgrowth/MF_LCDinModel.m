  
function [output]=MF_LCDinModel(ModelMatrix)
%Calculates the LCD and the Gradient of the LCD (is it stable enough?: just
%taking adjacent values) and stores it in the ModelMatrix.

%Test: ModelMatrix=CurrentModelMatrix;
   
% init LCD output measurement
matObjectCount=length(ModelMatrix.radius) ;  
cellMeasurement_Nuclei_LocalCellDensity = arrayfun(@(x) NaN(x,1),matObjectCount,'UniformOutput',false)';
%cellMeasurement_Nuclei_Edge = arrayfun(@(x) NaN(x,1),matObjectCount,'UniformOutput',false)';
%cellMeasurement_Nuclei_DistanceToEdge = arrayfun(@(x) NaN(x,1),matObjectCount,'UniformOutput',false)';
%cellMeasurement_Nuclei_Single = arrayfun(@(x) NaN(x,1),matObjectCount,'UniformOutput',false)';

%these are the values used to calculate the movie!
 intFilterSize = 486;
 intFilterSigma = intFilterSize*(25/150);
 intShrinkFactor = 13;
intImageBorderOverlap=8;
%round(sum(PSF(:))))=243


%%%%%%%%%%%%%%%%%%%PSF%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the point spread function (PSF) for image dilution
PSF = fspecial('gaussian',intFilterSize,intFilterSigma);
%adapt size of PSF (by factor 1/intShrinkFactor)
PSF = imresize(PSF,1/intShrinkFactor);

PSF = PSF - min(PSF(:));
PSF = PSF / max(PSF(:));
intMaxPSFValue = max(PSF(:));
intSingleCellFactor = 1.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
    % get all nuclei positions 
    matNucleiPositions = cat(1,ModelMatrix.Nuclei_Location);
       matNucleiPositions = ModelMatrix.Nuclei_Location;
    % shrink nuclei positions if necessary
    if intShrinkFactor~= 1
        matNucleiPositions = ceil(matNucleiPositions / intShrinkFactor);
    end
    % fix weird bug with nuclei positions of 0
    matNucleiPositions(matNucleiPositions==0)=1;
    matNucleiPositions(matNucleiPositions<0)=1;
        
    % create map with dots for each cell
    % perhaps work with dots, gaussian blurred...
    matImageMapWithDots = zeros(ceil(ModelMatrix.columnnumber / intShrinkFactor),ceil(ModelMatrix.rownumber / intShrinkFactor));
    for iCell = 1:size(matNucleiPositions,1)
        matImageMapWithDots(matNucleiPositions(iCell,2),matNucleiPositions(iCell,1)) = 1;
    end

    clear iCell
    
    matImageMap = matImageMapWithDots;
    matImageMap = imfilter(matImageMap,PSF,'symmetric','conv');

%Calculate Gradient    
[gradientx, gradienty]=gradient(matImageMap);


% store LCD, gradient per object, from image
ModelMatrix.LCD= arrayfun(@(x,y) matImageMap(x,y),matNucleiPositions(:,2),matNucleiPositions(:,1));
ModelMatrix.Gradientx= ...
            arrayfun(@(x,y) gradientx(x,y),matNucleiPositions(:,2),matNucleiPositions(:,1));
ModelMatrix.Gradienty= ...
            arrayfun(@(x,y) gradienty(x,y),matNucleiPositions(:,2),matNucleiPositions(:,1));
        
     

   clear matNucleiPositions PSF matObjectCount intFilterSize intFilterSigma intShrinkFactor
  clear matObjectCount
 

%Quality test: scatter(ModelMatrix.LCDnew, ModelMatrix.LCD)
output=ModelMatrix;
