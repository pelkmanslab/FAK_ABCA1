function [list] = MF_ListObjectsAroundObject(inputmatrix, nucleiconsidered, filterthreshold)
%determines which objects eg. nucleis are withhin the filterthreshold range
%of the centroid of the object "nucleiconsidered". For this the location of
%the other objects are saved in the inputmatrix.


%1. Filter based on size: get all "relevant nuclei" in surrounding
position_nucleicenter=inputmatrix.Nuclei_Location(nucleiconsidered,:);
lowerboundaryA =position_nucleicenter(1,1)-filterthreshold;
upperboundaryA =position_nucleicenter(1,1)+filterthreshold;
lowerboundaryB=position_nucleicenter(1,2)-filterthreshold;
upperboundaryB=position_nucleicenter(1,2)+filterthreshold;

list=find((inputmatrix.Nuclei_Location(:,1)>lowerboundaryA &inputmatrix.Nuclei_Location(:,1)<upperboundaryA) & (inputmatrix.Nuclei_Location(:,2)>lowerboundaryB &inputmatrix.Nuclei_Location(:,2)<upperboundaryB));
list=setdiff(list,nucleiconsidered);