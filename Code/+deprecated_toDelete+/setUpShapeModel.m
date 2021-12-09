%% 

% Authors:
%     Meghan Keast
%     Centre for Sport Research
%     Deakin University
%     mfkeast@deakin.edu.au
%
%     Aaron Fox
%     Centre for Sport Research
%     Deakin University
%     aaron.f@deakin.edu.au

% This function takes the segmented tibias and fibulas, aligns them to a
% global coordinate system, remeshes them to a consistent node number, and
% outputs them to a directory in the ShapeModel folder for use with the
% Python gias2 code. The steps included here essentially help out the shape
% model code work a little better and faster.

%%%%% TODO: consider extracting cortical:trabecular volume parameters as
%%%%% part of this function...

%% Set-up

%Set home directory
homeDir = cd;

%Add supplementary code path
addpath([pwd,'\Supplementary']);

%Navigate to segmentation directory
cd('..\Segmentation\');

%Grab the case names
f = dir();
for ff = 3:length(f)
    caseID{ff-2} = f(ff).name;
end

%Set options for remeshing
% % % optionStruct_tib.pointSpacing = 2.5; %Set desired point spacing -- another approach for remeshing
optionStruct_tib.nb_pts = 5000; %Set desired number of points
optionStruct_tib.disp_on = 0; % Turn off command window text display
% % % optionStruct_fib.pointSpacing = 2.5; %Set desired point spacing -- another approach for remeshing
optionStruct_fib.nb_pts = 2500; %Set desired number of points
optionStruct_fib.disp_on = 0; % Turn off command window text display

%Set options for visualising each remeshing
visualiseRemesh = false; %change to true to view each iteration

%Create space to store text strings for gias2 modelling
%Tibia
rbfregTib = {};
rigidregTib = {};
pcalistTib = {};
pcalistscaledTib = {};
%Trabecular
rbfregTrab = {};
rigidregTrab = {};
pcalistTrab = {};
pcalistscaledTrab = {};
%Tibia-fibula
rbfregTibFib = {};
rigidregTibFib = {};
pcalistTibFib = {};
pcalistscaledTibFib = {};

%% Extract and process surfaces

%Loop through cases
for caseInd = 1:length(caseID)
    
    %Navigate to case ID
    cd(caseID{caseInd});
    
    %Get case number
    caseNo = strsplit(caseID{caseInd},'-');
    caseNo = caseNo{2};
    
    %Create the check for no fibula in the case of case-134065
    if strcmp(caseNo, '134065')
        noFib = true;
    else
        noFib = false;
    end
    
    %Load the tibia - cortical
    %Uses GIBBON functionality
    [tibiaSTLstruct] = import_STL([caseNo,'-tibia-cortical-remesh.stl']);
    tibiaF = tibiaSTLstruct.solidFaces{1}; %Faces
    tibiaV = tibiaSTLstruct.solidVertices{1}; %Vertices
    [tibiaF,tibiaV] = mergeVertices(tibiaF,tibiaV);
    
    %Load the tibia - trabecular
    %Uses GIBBON functionality
    [trabSTLstruct] = import_STL([caseNo,'-tibia-trabecular.stl']);
    trabF = trabSTLstruct.solidFaces{1}; %Faces
    trabV = trabSTLstruct.solidVertices{1}; %Vertices
    [trabF,trabV] = mergeVertices(trabF,trabV);
    
    %Load the fibula
    if ~noFib
        %Uses GIBBON functionality
        [fibulaSTLstruct] = import_STL([caseNo,'-fibula.stl']);
        fibulaF = fibulaSTLstruct.solidFaces{1}; %Faces
        fibulaV = fibulaSTLstruct.solidVertices{1}; %Vertices
        [fibulaF,fibulaV] = mergeVertices(fibulaF,fibulaV);
    end

    %Try to remesh and continue with case if possible
    try
        
        %Run initial smoothing on surfaces
        smoothPar.n = 2; %2 iterations
        smoothPar.Method = 'LAP'; %Laplacian smoothing method        
        tibiaV = patchSmooth(tibiaF,tibiaV,[],smoothPar);
        trabV = patchSmooth(trabF,trabV,[],smoothPar);
        if ~noFib
            fibulaV = patchSmooth(fibulaF,fibulaV,[],smoothPar);
        end

        %Remesh using ggremesh to reduce complexity
        [tibiaF,tibiaV] = ggremesh(tibiaF, tibiaV, optionStruct_tib);
        
        %Check for matching point number and display if not
        if length(tibiaV) ~= optionStruct_tib.nb_pts
            fprintf('Number of tibia vertices does not match requested for case %s\n',caseNo);
        end
        
        %Remesh using ggremesh to reduce complexity
        [trabF,trabV] = ggremesh(trabF,trabV, optionStruct_tib);
        
        %Check for matching point number and display if not
        if length(trabV) ~= optionStruct_tib.nb_pts
            %First time here at least retry the remesh again
            [trabF,trabV] = ggremesh(trabF,trabV, optionStruct_tib);
            %Check again
            if length(trabV) ~= optionStruct_tib.nb_pts
                %Print out summary this time
                fprintf('Number of trabecular vertices does not match requested for case %s\n',caseNo);
            end
        end
        
        %Remesh using ggremesh to reduce complexity
        if ~noFib
            [fibulaF,fibulaV] = ggremesh(fibulaF,fibulaV, optionStruct_fib);
        
            %Check for matching point number and display if not
            if length(fibulaV) ~= optionStruct_fib.nb_pts
                %First time here at least retry the remesh again
                [fibulaF,fibulaV] = ggremesh(fibulaF,fibulaV, optionStruct_fib);
                %Check again
                if length(fibulaV) ~= optionStruct_fib.nb_pts
                    %Print out summary this time
                    fprintf('Number of fibula vertices does not match requested for case %s\n',caseNo);
                end
            end
        end
        
        %Check for holes in surfaces
        %This shouldn't occur with the remeshing, but check for safety...
        if ~isempty(patchBoundary(tibiaF, tibiaV))
            fprintf('Holes in tibia mesh for %s\n...',caseNo)
        end
        if ~isempty(patchBoundary(trabF, trabV))
            fprintf('Holes in trabecular mesh for %s\n...',caseNo)
        end
        if ~noFib
            if ~isempty(patchBoundary(fibulaF, fibulaV))
                fprintf('Holes in fibula mesh for %s\n...',caseNo)
            end
        end
               
        %Import tibial landmarks
        tibPoints = [{'LC'},{'LM'},{'MC'},{'MM'}];
        for pp = 1:length(tibPoints)
            tree = xml_read([tibPoints{pp},'.txt']);
            landmarks.(char(tree.Point.Name)) = tree.Point.Coordinate;% / 1000;
            clear tree 
        end
        clear pp

        %Create new landmarks
        landmarks.IM = midPoint3d(landmarks.LM,landmarks.MM);
        landmarks.IC = midPoint3d(landmarks.LC,landmarks.MC);

        %Create and align tibial coordinate system
        if ~noFib
            [tibiaF, tibiaV, fibulaF, fibulaV, trabF, trabV] = ...
                rotateSurfacesGlobal(landmarks, noFib, tibiaF, tibiaV, fibulaF, fibulaV, trabF, trabV);
        elseif noFib
            [tibiaF, tibiaV, ~, ~, trabF, trabV] = ...
                rotateSurfacesGlobal(landmarks, noFib, tibiaF, tibiaV, [], [], trabF, trabV);
        end

        %Visualise
        if visualiseRemesh
            %Create figure
            h = cFigure; hold on
            %Plot surfaces
            gpatch(tibiaF,tibiaV,'gw','k', 0.3); %low alpha to see trabecular 
            gpatch(trabF,trabV,'rw','k');
            if ~noFib
                gpatch(fibulaF,fibulaV,'bw','k');
            end
            axisGeom; camlight headlight
            title('Imported and Aligned Tibia-Fibula');
            %Pause
            pause
            %Close figure
            close all
        end

        %Export the STL
               
        %Merge objects into one surface
        if ~noFib
            mergedSurfaceV = [tibiaV; fibulaV];
            mergedSurfaceF = [tibiaF; fibulaF+length(tibiaV)];
        end

        %Create the stl structures
        %Tibia
        stlStructTib.solidNames = {['tibia-',caseNo]}; %names of parts
        stlStructTib.solidVertices = {tibiaV}; %Vertices
        stlStructTib.solidFaces = {tibiaF}; %Faces
        stlStructTib.solidNormals={[]};
        %Tibia-Fibula
        if ~noFib
            stlStructTibFib.solidNames = {['tibia-fibula-',caseNo]}; %names of parts
            stlStructTibFib.solidVertices = {mergedSurfaceV}; %Vertices
            stlStructTibFib.solidFaces = {mergedSurfaceF}; %Faces
            stlStructTibFib.solidNormals={[]};
        end
        %Trabecular
        stlStructTrab.solidNames = {['trabecular-',caseNo]}; %names of parts
        stlStructTrab.solidVertices = {trabV}; %Vertices
        stlStructTrab.solidFaces = {trabF}; %Faces
        stlStructTrab.solidNormals={[]};

        %Export the stls
        export_STL_txt(['..\..\ShapeModels\tibia\segmentations\',caseNo,'-tibia.stl'], stlStructTib);
        if ~noFib
            export_STL_txt(['..\..\ShapeModels\tibia-fibula\segmentations\',caseNo,'-tibia-fibula.stl'], stlStructTibFib);
        end
        export_STL_txt(['..\..\ShapeModels\trabecular\segmentations\',caseNo,'-trabecular.stl'], stlStructTrab);
        
        %Write the case details to a string array for the pca text files
        %Tibia
        rbfregTib{length(rbfregTib)+1,1} = ['segmentations/',caseNo,'-tibia.stl'];
        rigidregTib{length(rigidregTib)+1,1} = ['fitted_meshes/',caseNo,'-tibia_rbfreg.stl'];
        pcalistTib{length(pcalistTib)+1,1} = ['aligned_meshes/',caseNo,'-tibia_rbfreg_rigidreg.stl'];
		pcalistscaledTib{length(pcalistscaledTib)+1,1} = ['aligned_meshes_scaled/',caseNo,'-tibia_rbfreg_rigidreg.stl'];
        %Tibia-Fibula
        if ~noFib
            rbfregTibFib{length(rbfregTibFib)+1,1} = ['segmentations/',caseNo,'-tibia-fibula.stl'];
            rigidregTibFib{length(rigidregTibFib)+1,1} = ['fitted_meshes/',caseNo,'-tibia-fibula_rbfreg.stl'];
            pcalistTibFib{length(pcalistTibFib)+1,1} = ['aligned_meshes/',caseNo,'-tibia-fibula_rbfreg_rigidreg.stl'];
			pcalistscaledTibFib{length(pcalistscaledTibFib)+1,1} = ['aligned_meshes_scaled/',caseNo,'-tibia-fibula_rbfreg_rigidreg.stl'];
        end
        %Trabecular
        rbfregTrab{length(rbfregTrab)+1,1} = ['segmentations/',caseNo,'-trabecular.stl'];
        rigidregTrab{length(rigidregTrab)+1,1} = ['fitted_meshes/',caseNo,'-trabecular_rbfreg.stl'];
        pcalistTrab{length(pcalistTrab)+1,1} = ['aligned_meshes/',caseNo,'-trabecular_rbfreg_rigidreg.stl'];
		pcalistscaledTrab{length(pcalistscaledTrab)+1,1} = ['aligned_meshes_scaled/',caseNo,'-trabecular_rbfreg_rigidreg.stl'];
        
        %Determine cortical:trabecular volume along bone
        
        %%%%% TODO: SHELVING THIS FOR HERE - AS I DON'T THINK VOLUME EVERY
        %%%%% 5% IS APPROPIRATE, IT'S MORE LIKE THE THICKNESS AT EACH SORT
        %%%%% OF STEP UP THE TIBIA IS APPROPRIATE - I.E. CUT THE TIBIA
        %%%%% EVERY 1% AND ASSESS THE CORTICAL THICKNESS AT THAT CUT
        %%%%% LEVEL...?
        
% % %         calcRelativeTrabecularVol(tibiaF, tibiaV, trabF, trabV);            

        %Disp progress
        disp(['Case ',caseNo,' successfully remeshed and aligned.']);
        
    catch
        
        %Display catch issue
        disp(['Unable to remesh case ',caseNo,'. Skipped...'])
    
    end
    
    %Clear variables for re-loop
    clearvars -except visualiseRemesh caseID homeDir caseInd optionStruct_tib optionStruct_fib pcalistTib pcalistscaledTib rbfregTib rigidregTib pcalistTrab pcalistscaledTrab rbfregTrab rigidregTrab pcalistTibFib pcalistscaledTibFib rbfregTibFib rigidregTibFib
    
    %Navigate back to segmentation directory
    cd('..');
    
end

%Before writing these case details to text files, we set the top row to
%correspond to a mesh we want all others registered against (i.e. a target
%mesh). For this we use a case that is closest to the mean surface in line
%with existing work (i.e. Bruce et al. 2021, Computer Methods Biomech
%Biomed Eng, doi: 10.1080/10255842.2021.1985111). Some testing revealed
%that the best use case for this was case-147211, hence we rewrite this
%case to sit at the top of all the file lists. 

%Set reference case number
refCase = '147211';

%Reorder the arrays
%Tibia
rbfregTib = [rbfregTib(find(contains(rbfregTib, refCase)));
    rbfregTib(1:find(contains(rbfregTib, refCase))-1);
    rbfregTib(find(contains(rbfregTib, refCase))+1:end)];
rigidregTib = [rigidregTib(find(contains(rigidregTib, refCase)));
    rigidregTib(1:find(contains(rigidregTib, refCase))-1);
    rigidregTib(find(contains(rigidregTib, refCase))+1:end)];
pcalistTib = [pcalistTib(find(contains(pcalistTib, refCase)));
    pcalistTib(1:find(contains(pcalistTib, refCase))-1);
    pcalistTib(find(contains(pcalistTib, refCase))+1:end)];
pcalistscaledTib = [pcalistscaledTib(find(contains(pcalistscaledTib, refCase)));
    pcalistscaledTib(1:find(contains(pcalistscaledTib, refCase))-1);
    pcalistscaledTib(find(contains(pcalistscaledTib, refCase))+1:end)];
%Trabecular
rbfregTrab = [rbfregTrab(find(contains(rbfregTrab, refCase)));
    rbfregTrab(1:find(contains(rbfregTrab, refCase))-1);
    rbfregTrab(find(contains(rbfregTrab, refCase))+1:end)];
rigidregTrab = [rigidregTrab(find(contains(rigidregTrab, refCase)));
    rigidregTrab(1:find(contains(rigidregTrab, refCase))-1);
    rigidregTrab(find(contains(rigidregTrab, refCase))+1:end)];
pcalistTrab = [pcalistTrab(find(contains(pcalistTrab, refCase)));
    pcalistTrab(1:find(contains(pcalistTrab, refCase))-1);
    pcalistTrab(find(contains(pcalistTrab, refCase))+1:end)];
pcalistscaledTrab = [pcalistscaledTrab(find(contains(pcalistscaledTrab, refCase)));
    pcalistscaledTrab(1:find(contains(pcalistscaledTrab, refCase))-1);
    pcalistscaledTrab(find(contains(pcalistscaledTrab, refCase))+1:end)];
%Tibia-fibula
rbfregTibFib = [rbfregTibFib(find(contains(rbfregTibFib, refCase)));
    rbfregTibFib(1:find(contains(rbfregTibFib, refCase))-1);
    rbfregTibFib(find(contains(rbfregTibFib, refCase))+1:end)];
rigidregTibFib = [rigidregTibFib(find(contains(rigidregTibFib, refCase)));
    rigidregTibFib(1:find(contains(rigidregTibFib, refCase))-1);
    rigidregTibFib(find(contains(rigidregTibFib, refCase))+1:end)];
pcalistTibFib = [pcalistTibFib(find(contains(pcalistTibFib, refCase)));
    pcalistTibFib(1:find(contains(pcalistTibFib, refCase))-1);
    pcalistTibFib(find(contains(pcalistTibFib, refCase))+1:end)];
pcalistscaledTibFib = [pcalistscaledTibFib(find(contains(pcalistscaledTibFib, refCase)));
    pcalistscaledTibFib(1:find(contains(pcalistscaledTibFib, refCase))-1);
    pcalistscaledTibFib(find(contains(pcalistscaledTibFib, refCase))+1:end)];

%Write the case details to the text files needed for the shape model
writeShapeModelTextFiles(rbfregTib, rigidregTib, pcalistTib, pcalistscaledTib, ...
    rbfregTibFib, rigidregTibFib, pcalistTibFib, pcalistscaledTibFib, ...
    rbfregTrab, rigidregTrab, pcalistTrab, pcalistscaledTrab);

%% Finish up

%Return to home directory
cd(homeDir);

%% ----- End of shapeModelSetUp.m ----- %%