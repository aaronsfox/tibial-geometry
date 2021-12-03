function writeShapeModelTextFiles(rbfregTib, rigidregTib, pcalistTib, pcalistscaledTib, ...
    rbfregTibFib, rigidregTibFib, pcalistTibFib, pcalistscaledTibFib, ...
    rbfregTrab, rigidregTrab, pcalistTrab, pcalistscaledTrab)

    %% Function that cleans up set-up script by shifting code to write the
    %  text files for the shape model functions to this separate script.
    
    %Start by removing any old list files
    %Tibia
    %Remove rbfreg file if present
    if isfile('..\ShapeModels\tibia\rbfreg_list.txt')
        delete('..\ShapeModels\tibia\rbfreg_list.txt')
    end
    %Remove rigidreg if present
    if isfile('..\ShapeModels\tibia\rigidreg_list.txt')
        delete('..\ShapeModels\tibia\rigidreg_list.txt')
    end
    %Remove pca list if present
    if isfile('..\ShapeModels\tibia\pca_list.txt')
        delete('..\ShapeModels\tibia\pca_list.txt')
    end
    %Remove pca scaled list if present
    if isfile('..\ShapeModels\tibia\pca_list_scaled.txt')
        delete('..\ShapeModels\tibia\pca_list_scaled.txt')
    end
    
    %Tibia-Fibula
    %Remove rbfreg file if present
    if isfile('..\ShapeModels\tibia-fibula\rbfreg_list.txt')
        delete('..\ShapeModels\tibia-fibula\rbfreg_list.txt')
    end
    %Remove rigidreg if present
    if isfile('..\ShapeModels\tibia-fibula\rigidreg_list.txt')
        delete('..\ShapeModels\tibia-fibula\rigidreg_list.txt')
    end
    %Remove pca list if present
    if isfile('..\ShapeModels\tibia-fibula\pca_list.txt')
        delete('..\ShapeModels\tibia-fibula\pca_list.txt')
    end
    %Remove pca scaled list if present
    if isfile('..\ShapeModels\tibia-fibula\pca_list_scaled.txt')
        delete('..\ShapeModels\tibia-fibula\pca_list_scaled.txt')
    end
    
    %Trabecular
    %Remove rbfreg file if present
    if isfile('..\ShapeModels\trabecular\rbfreg_list.txt')
        delete('..\ShapeModels\trabecular\rbfreg_list.txt')
    end
    %Remove rigidreg if present
    if isfile('..\ShapeModels\trabecular\rigidreg_list.txt')
        delete('..\ShapeModels\trabecular\rigidreg_list.txt')
    end
    %Remove pca list if present
    if isfile('..\ShapeModels\trabecular\pca_list.txt')
        delete('..\ShapeModels\trabecular\pca_list.txt')
    end
    %Remove pca scaled list if present
    if isfile('..\ShapeModels\trabecular\pca_list_scaled.txt')
        delete('..\ShapeModels\trabecular\pca_list_scaled.txt')
    end

    %Write the created list sets to file
    %Tibia
    %rbfreg
    fid = fopen('..\ShapeModels\tibia\rbfreg_list.txt','wt');
    for caseNo = 1:length(rbfregTib)
        %Check if current case was successful
        if ~isempty(rbfregTib{caseNo,1})
            %Write to file
            fprintf(fid, rbfregTib{caseNo,1});
            %Check to see if new line is needed
            if caseNo < length(rbfregTib)
                fprintf(fid,'\n');
            end
        end
    end
    fclose(fid);
    %rigidreg
    fid = fopen('..\ShapeModels\tibia\rigidreg_list.txt','wt');
    for caseNo = 1:length(rigidregTib)
        %Check if current case was successful
        if ~isempty(rigidregTib{caseNo,1})
            %Write to file
            fprintf(fid, rigidregTib{caseNo,1});
            %Check to see if new line is needed
            if caseNo < length(rigidregTib)
                fprintf(fid,'\n');
            end
        end
    end
    fclose(fid);
    %pca
    fid = fopen('..\ShapeModels\tibia\pca_list.txt','wt');
    for caseNo = 1:length(pcalistTib)
        %Check if current case was successful
        if ~isempty(pcalistTib{caseNo,1})
            %Write to file
            fprintf(fid, pcalistTib{caseNo,1});
            %Check to see if new line is needed
            if caseNo < length(pcalistTib)
                fprintf(fid,'\n');
            end
        end
    end
    fclose(fid);
    %pca scaled
    fid = fopen('..\ShapeModels\tibia\pca_list_scaled.txt','wt');
    for caseNo = 1:length(pcalistscaledTib)
        %Check if current case was successful
        if ~isempty(pcalistscaledTib{caseNo,1})
            %Write to file
            fprintf(fid, pcalistscaledTib{caseNo,1});
            %Check to see if new line is needed
            if caseNo < length(pcalistscaledTib)
                fprintf(fid,'\n');
            end
        end
    end
    fclose(fid);

    %Tibia-Fibula
    %rbfreg
    fid = fopen('..\ShapeModels\tibia-fibula\rbfreg_list.txt','wt');
    for caseNo = 1:length(rbfregTibFib)
        %Check if current case was successful
        if ~isempty(rbfregTibFib{caseNo,1})
            %Write to file
            fprintf(fid, rbfregTibFib{caseNo,1});
            %Check to see if new line is needed
            if caseNo < length(rbfregTibFib)
                fprintf(fid,'\n');
            end
        end
    end
    fclose(fid);
    %rigidreg
    fid = fopen('..\ShapeModels\tibia-fibula\rigidreg_list.txt','wt');
    for caseNo = 1:length(rigidregTibFib)
        %Check if current case was successful
        if ~isempty(rigidregTibFib{caseNo,1})
            %Write to file
            fprintf(fid, rigidregTibFib{caseNo,1});
            %Check to see if new line is needed
            if caseNo < length(rigidregTibFib)
                fprintf(fid,'\n');
            end
        end
    end
    fclose(fid);
    %pca
    fid = fopen('..\ShapeModels\tibia-fibula\pca_list.txt','wt');
    for caseNo = 1:length(pcalistTibFib)
        %Check if current case was successful
        if ~isempty(pcalistTibFib{caseNo,1})
            %Write to file
            fprintf(fid, pcalistTibFib{caseNo,1});
            %Check to see if new line is needed
            if caseNo < length(pcalistTibFib)
                fprintf(fid,'\n');
            end
        end
    end
    fclose(fid);
    %pca scaled
    fid = fopen('..\ShapeModels\tibia-fibula\pca_list_scaled.txt','wt');
    for caseNo = 1:length(pcalistscaledTibFib)
        %Check if current case was successful
        if ~isempty(pcalistscaledTibFib{caseNo,1})
            %Write to file
            fprintf(fid, pcalistscaledTibFib{caseNo,1});
            %Check to see if new line is needed
            if caseNo < length(pcalistscaledTibFib)
                fprintf(fid,'\n');
            end
        end
    end
    fclose(fid);
    
    %Trabecular
    %rbfreg
    fid = fopen('..\ShapeModels\trabecular\rbfreg_list.txt','wt');
    for caseNo = 1:length(rbfregTrab)
        %Check if current case was successful
        if ~isempty(rbfregTrab{caseNo,1})
            %Write to file
            fprintf(fid, rbfregTrab{caseNo,1});
            %Check to see if new line is needed
            if caseNo < length(rbfregTrab)
                fprintf(fid,'\n');
            end
        end
    end
    fclose(fid);
    %rigidreg
    fid = fopen('..\ShapeModels\trabecular\rigidreg_list.txt','wt');
    for caseNo = 1:length(rigidregTrab)
        %Check if current case was successful
        if ~isempty(rigidregTrab{caseNo,1})
            %Write to file
            fprintf(fid, rigidregTrab{caseNo,1});
            %Check to see if new line is needed
            if caseNo < length(rigidregTrab)
                fprintf(fid,'\n');
            end
        end
    end
    fclose(fid);
    %pca
    fid = fopen('..\ShapeModels\trabecular\pca_list.txt','wt');
    for caseNo = 1:length(pcalistTrab)
        %Check if current case was successful
        if ~isempty(pcalistTrab{caseNo,1})
            %Write to file
            fprintf(fid, pcalistTrab{caseNo,1});
            %Check to see if new line is needed
            if caseNo < length(pcalistTrab)
                fprintf(fid,'\n');
            end
        end
    end
    fclose(fid);
    %pca scaled
    fid = fopen('..\ShapeModels\trabecular\pca_list_scaled.txt','wt');
    for caseNo = 1:length(pcalistscaledTrab)
        %Check if current case was successful
        if ~isempty(pcalistscaledTrab{caseNo,1})
            %Write to file
            fprintf(fid, pcalistscaledTrab{caseNo,1});
            %Check to see if new line is needed
            if caseNo < length(pcalistscaledTrab)
                fprintf(fid,'\n');
            end
        end
    end
    fclose(fid);

end