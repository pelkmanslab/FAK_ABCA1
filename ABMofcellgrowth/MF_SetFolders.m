function[foldername1,foldername2,foldername3,foldername4,foldername5]=MF_SetFolders(foldername)


foldername1=[foldername '\\images'];
foldername2=[foldername '\\CurrentModelMatrices'];
foldername3=[foldername '\\interactionMatrices'];
foldername4=[foldername '\\energymatrices'];
foldername5=[foldername '\\images_2'];

% create output directory if it doesn't exist.
    if ~fileattrib(foldername)
        mkdir(foldername)
    end
    
    if ~fileattrib(foldername1)
       mkdir(foldername1)
    end

    if ~fileattrib(foldername2)
        mkdir(foldername2)
    end

    if ~fileattrib(foldername3)
        mkdir(foldername3)
    end

    if ~fileattrib(foldername4)
        mkdir(foldername4)
    end
    
    if ~fileattrib(foldername5)
        mkdir(foldername5)
    end

   

