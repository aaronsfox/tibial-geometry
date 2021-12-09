

   




%% Set-up

%Set home directory
homeDir = cd;

%Turn warnings off
warning off

%Add supplementary code path
addpath(genpath([pwd,'\Supplementary']));

%Navigate to segmentation directory
cd('..\Segmentation\');

%Grab the case names
f = dir();
for ff = 3:length(f)
    caseID{ff-2} = f(ff).name;
end

%Set options for remeshing
% % % optionStruct_tib.pointSpacing = 2.5; %Set desired point spacing -- another approach for remeshing
% % % optionStruct_tib.nb_pts = 5000; %Set desired number of points
optionStruct_tib.nb_pts = 3500; %Set desired number of points
optionStruct_tib.disp_on = 0; % Turn off command window text display
% % % optionStruct_fib.pointSpacing = 2.5; %Set desired point spacing -- another approach for remeshing
% % % optionStruct_fib.nb_pts = 2500; %Set desired number of points
optionStruct_fib.nb_pts = 2000; %Set desired number of points
optionStruct_fib.disp_on = 0; % Turn off command window text display

%Set options for CPD algorithm
optCPD.method = 'nonrigid'; % use nonrigid registration
optCPD.beta = 2;            % the width of Gaussian kernel (smoothness)
optCPD.lambda = 3;          % regularization weight
optCPD.viz = 0;             % don't visualise
optCPD.outliers = 0;        % don't account for outliers
optCPD.fgt = 0;             % do not use FGT (default)
optCPD.normalize = 1;       % normalize to unit variance and zero mean before registering (default)
optCPD.corresp = 0;         % compute correspondence vector at the end of registration (not being estimated by default)
optCPD.max_it = 100;        % max number of iterations
optCPD.tol = 1e-4;          % tolerance
optCPD.corresp = 0;         % estimate correspondence

%Set options for visualising each remeshing
visualiseRemesh = false; %change to true to view each iteration

%Set threshold for keeping PCs of models (i.e. % variance explained)
varExpThreshold = 95;

%% Remesh and align surfaces

%All surface meshes need to be aligned to a single target surface. For this
%we use a case that is closest to the mean surface in line with existing
% work (i.e. Bruce et al. 2021, Computer Methods Biomech Biomed Eng, doi:
%10.1080/10255842.2021.1985111). Some testing revealed that the best use
%case for this was case-147211, hence we ensure this is processed first and
%all subsequent meshes are aligned to this

%Set reference case number
refCase = '147211';

%Reorder case list to make target mesh first
caseID = [caseID(find(contains(caseID, refCase))), ...
    caseID(1:find(contains(caseID, refCase))-1), ...
    caseID(find(contains(caseID, refCase))+1:end)];

%Remesh and align all cases

%Create waitbar to monitor progress
wbar = waitbar(0,'Remeshing and aligning surfaces....');

%Loop through cases
for caseInd = 1:length(caseID)
    
    %Navigate to case ID
    cd(caseID{caseInd});
    
    %Get case number
    caseNo = strsplit(caseID{caseInd},'-');
    caseNo = caseNo{2};
    
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
    %Uses GIBBON functionality
    [fibulaSTLstruct] = import_STL([caseNo,'-fibula.stl']);
    fibulaF = fibulaSTLstruct.solidFaces{1}; %Faces
    fibulaV = fibulaSTLstruct.solidVertices{1}; %Vertices
    [fibulaF,fibulaV] = mergeVertices(fibulaF,fibulaV);

    %Remesh surfaces
        
    %Run initial smoothing on surfaces
    smoothPar.n = 2; %2 iterations
    smoothPar.Method = 'LAP'; %Laplacian smoothing method        
    tibiaV = patchSmooth(tibiaF,tibiaV,[],smoothPar);
    trabV = patchSmooth(trabF,trabV,[],smoothPar);
    fibulaV = patchSmooth(fibulaF,fibulaV,[],smoothPar);

    %Remesh using ggremesh to reduce complexity
    [tibiaF,tibiaV] = ggremesh(tibiaF, tibiaV, optionStruct_tib);

    %Check for matching point number and display if not
    if length(tibiaV) ~= optionStruct_tib.nb_pts
        error('Number of tibia vertices does not match requested for case %s\n',caseNo);
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
            error('Number of trabecular vertices does not match requested for case %s\n',caseNo);
        end
    end

    %Remesh using ggremesh to reduce complexity
    [fibulaF,fibulaV] = ggremesh(fibulaF,fibulaV, optionStruct_fib);

    %Check for matching point number and display if not
    if length(fibulaV) ~= optionStruct_fib.nb_pts
        %First time here at least retry the remesh again
        [fibulaF,fibulaV] = ggremesh(fibulaF,fibulaV, optionStruct_fib);
        %Check again
        if length(fibulaV) ~= optionStruct_fib.nb_pts
            %Print out summary this time
            error('Number of fibula vertices does not match requested for case %s\n',caseNo);
        end
    end

    %Check for holes in surfaces
    %This shouldn't occur with the remeshing, but check for safety...
    if ~isempty(patchBoundary(tibiaF, tibiaV))
        error('Holes in tibia mesh for %s\n...',caseNo)
    end
    if ~isempty(patchBoundary(trabF, trabV))
        error('Holes in trabecular mesh for %s\n...',caseNo)
    end
    if ~isempty(patchBoundary(fibulaF, fibulaV))
        error('Holes in fibula mesh for %s\n...',caseNo)
    end

    %Import tibial landmarks
    %%%%% TODO: convert these xml into a csv for ease of use load in...
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
    [tibiaF, tibiaV, fibulaF, fibulaV, trabF, trabV] = ...
        rotateSurfacesGlobal(landmarks, tibiaF, tibiaV, fibulaF, fibulaV, trabF, trabV);

    %Visualise
    if visualiseRemesh
        %Create figure
        h = cFigure; hold on
        %Plot surfaces
        gpatch(tibiaF,tibiaV,'gw','k', 0.3); %low alpha to see trabecular 
        gpatch(trabF,trabV,'rw','k');
        gpatch(fibulaF,fibulaV,'bw','k');
        axisGeom; camlight headlight
        title('Imported and Aligned Tibia-Fibula');
        %Pause
        pause
        %Close figure
        close all
    end

    %If not the target mesh, align the surfaces
    if caseInd > 1
        
        %Tibia
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['APPLYING CPD ALGORITHM TO TIBIA OF CASE-',caseNo,'...']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Perform non-rigid registration using CPD to align surfaces
        [cpdTformTib] = cpd_register(data.tibia.V(:,:,1), tibiaV, optCPD);
        
        
        %%%%% CORRESPONDING WORKS SLIGHTLY DIFFERENTLY BUT PRETTY MUCH SAME
        %%%%% RESULT, STILL PICKS UP SHARED POINTS...
        
        %%% Get nearest 10 or so points, sort by distance and allocate
        %%% those points, then take out the point that's been used already?
        
        
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %         disp(['APPLYING ICP ALGORITHM TO TIBIA OF CASE-',caseNo,'...']);
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %         
% % %         %%%% TODO: could increase to 20 iterations?????
% % %         
% % %         %Perform non-rigid registration using the ICP algorithm to assist in
% % %         %identifying corresponding points.
% % % % % %         [registeredPtsTib] = nonrigidICPv1(data.tibia.V(:,:,1), movingRegTib_cpd.Location, ...
% % % % % %             data.tibia.F(:,:,1), tibiaF, 10, 1);
% % %         [registeredPtsTib] = nonrigidICPv1(data.tibia.V(:,:,1), tibiaV, ...
% % %             data.tibia.F(:,:,1), tibiaF, 10, 0);

        %%%%% TODO: what about knnsearch finding the same points as nearest
        %%%%% neighbours? need to include distance and find closest for
        %%%%% each point without including same? Reducing mesh density
        %%%%% might help with this...?

        %Identify the matching points in the target mesh against the
        %registration to identify the corresponding point indices
        regSortIdxTib = knnsearch(cpdTformTib.Y, data.tibia.V(:,:,1));
        
        %Sort the registered points so they align with the reference mesh

        %%%% This seems appropriate --- but need to consider that fibula
        %%%% and trabecular are being rigidly aligned according to the
        %%%% tibia --- appropriate???
        
% % %         tibiaV_sorted = movingRegTib_cpd.Location(regSortIdxTib,:);
        tibiaV_sorted = tibiaV(regSortIdxTib,:);
        
        %Visualise to confirm that points now correspond. This can be done
        %by using the target mesh faces against the registered moving
        %points. The surface meshes should essentially come out looking the
        %same (i.e. overlapped) despite using different faces
        cFigure; hold on;
        subplot(1,3,1); hold on;
% % %         gpatch(data.tibia.F(:,:,1), data.tibia.V(:,:,1), 'gw', 'k', 0.3);
        gpatch(tibiaF, tibiaV, 'rw', 'none', 0.3);
        gpatch(data.tibia.F(:,:,1), tibiaV_sorted, 'gw', 'k');
        title('Registered Tibia');
        axisGeom;
        
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %         disp(['APPLYING PROCRUSTES TO TIBIA OF CASE-',caseNo,'...']);
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %         
% % %         %Perform rigid transformation using procrustes analysis
% % %         [tformTib_proc, movingRegTib_proc] = procrustes(data.tibia.V(:,:,1), ...
% % %             tibiaV_sorted, 'scaling', false);
        
        %Allocate to original variables
        tibiaF = data.tibia.F(:,:,1); %use target mesh faces
        tibiaV = tibiaV_sorted; %use sorted registered points
        
        %Trabecular
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['APPLYING CPD ALGORITHM TO TRABECULAR OF CASE-',caseNo,'...']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Perform non-rigid registration using CPD to align surfaces
        [cpdTformTrab] = cpd_register(data.trab_tibia.V(:,:,1), trabV, optCPD);
        
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %         disp(['APPLYING ICP ALGORITHM TO TRABECULAR OF CASE-',caseNo,'...']);
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %         
% % %         %Perform non-rigid registration using the ICP algorithm to assist in
% % %         %identifying corresponding points.
% % %         [registeredPtsTrab] = nonrigidICPv1(data.trab.V(:,:,1), trabV, ...
% % %             data.trab.F(:,:,1), trabF, 10, 0);

        %Identify the matching points in the target mesh against the
        %registration to identify the corresponding point indices
        regSortIdxTrab = knnsearch(cpdTformTrab.Y, data.trab_tibia.V(:,:,1));
        
        %Sort the registered points so they align with the shape model mean points
        trabV_sorted = trabV(regSortIdxTrab,:);
        
        %Visualise to confirm that points now correspond. This can be done
        %by using the target mesh faces against the registered moving
        %points. The surface meshes should essentially come out looking the
        %same (i.e. overlapped) despite using different faces
% % %         cFigure; hold on
        subplot(1,3,2); hold on;
% % %         gpatch(data.trab.F(:,:,1), data.trab.V(:,:,1), 'gw', 'k', 0.3);
        gpatch(trabF, trabV, 'rw', 'none', 0.3);
        gpatch(data.trab_tibia.F(:,:,1), trabV_sorted, 'gw', 'k');
        title('Registered Trabecular');
        axisGeom;
        
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %         disp(['APPLYING PROCRUSTES TO TRABECULAR OF CASE-',caseNo,'...']);
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %         
% % %         %Perform rigid transformation using procrustes analysis
% % %         [tformTrab_proc, movingRegTrab_proc] = procrustes(data.trab.V(:,:,1), ...
% % %             trabV_sorted, 'scaling', false);
        
        %Allocate to original variables
        trabF = data.trab_tibia.F(:,:,1); %use target mesh faces
        trabV = trabV_sorted; %use sorted registered points
        
        %Fibula

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['APPLYING CPD ALGORITHM TO FIBULA OF CASE-',caseNo,'...']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%% UP TO HERE SHOULD WORK FOR NOW....

        %Perform non-rigid registration using CPD to align surfaces
        [cpdTformFib] = cpd_register(data.fibula.V(:,:,1), fibulaV, optCPD);

% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %         disp(['APPLYING ICP ALGORITHM TO FIBULA OF CASE-',caseNo,'...']);
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % %         %Perform non-rigid registration using the ICP algorithm to assist in
% % %         %identifying corresponding points.
% % %         [registeredPtsFib] = nonrigidICPv1(data.fibula.V(:,:,1), fibulaV, ...
% % %             data.fibula.F(:,:,1), fibulaF, 10, 0);

        %Identify the matching points in the target mesh against the
        %registration to identify the corresponding point indices
        regSortIdxFib = knnsearch(cpdTformFib.Y, data.fibula.V(:,:,1));

        %Sort the registered points so they align with the shape model mean points
        fibulaV_sorted = fibulaV(regSortIdxFib,:);

        %Visualise to confirm that points now correspond. This can be done
        %by using the target mesh faces against the registered moving
        %points. The surface meshes should essentially come out looking the
        %same (i.e. overlapped) despite using different faces
% % %         cFigure; hold on
        subplot(1,3,3); hold on;
% % %         gpatch(data.fibula.F(:,:,1), data.fibula.V(:,:,1), 'gw', 'k', 0.3);
        gpatch(fibulaF, fibulaV, 'rw', 'none', 0.3);
        gpatch(data.fibula.F(:,:,1), fibulaV_sorted, 'gw', 'k');
        title('Registered Fibula');
        axisGeom;

% % %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %             disp(['APPLYING PROCRUSTES TO TIBIA OF CASE-',caseNo,'...']);
% % %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % %             %Perform rigid transformation using procrustes analysis
% % %             [tformFib_proc, movingRegFib_proc] = procrustes(data.fibula.V(:,:,1), ...
% % %                 fibulaV_sorted, 'scaling', false);

        %Allocate to original variables
        fibulaF = data.fibula.F(:,:,1); %use target mesh faces
        fibulaV = fibulaV_sorted; %use sorted registered points

        %Merge the tibia and fibula vertices for rigid alignment
        tibia_fibulaV = [tibiaV; fibulaV];
        tibia_fibulaF = data.tibia_fibula.F(:,:,1);
        
        %Export the subplot figure
        export_fig([caseNo,'_registeredSurfaces.png'],'-m1');
        close

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['APPLYING PROCRUSTES TO SURFACES OF CASE-',caseNo,'...']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Align based on tibia only first

        %Perform rigid transformation using procrustes analysis
        [~, regTibV, tformTib] = procrustes(data.tibia.V(:,:,1), ...
            tibiaV, 'scaling', false);
        
        %Visualise transform
        cFigure; hold on;
        plotV(data.tibia.V(:,:,1), 'k.');
        plotV(tibiaV, 'r.');
        plotV(regTibV, 'g.');
        axisGeom;
        legend('Reference','Original','Transformed');
        export_fig([caseNo,'_rigidTransformTibia.png'],'-m1');
        close
        
        %Replace original with transformed points
        tibiaV = regTibV;
        
        %Apply same transform to trabecular
        trab_tibiaV = trabV * tformTib.T + mean(tformTib.c);
        trab_tibiaF = trabF;
        
        %Align based on tibia and fibula second
        
        %Perform rigid transformation using procrustes analysis
        [~, regTibFibV, tformTibFib] = procrustes(data.tibia_fibula.V(:,:,1), ...
            tibia_fibulaV, 'scaling', false);
        
        %Visualise transform
        cFigure; hold on;
        plotV(data.tibia_fibula.V(:,:,1), 'k.');
        plotV(tibia_fibulaV, 'r.');
        plotV(regTibFibV, 'g.');
        axisGeom;
        legend('Reference','Original','Transformed');
        export_fig([caseNo,'_rigidTransformTibiaFibula.png'],'-m1');
        close
        
        %Replace original with transformed points
        tibia_fibulaV = regTibFibV;
        
        %Apply same transform to trabecular
        trab_tibia_fibulaV = trabV * tformTibFib.T + mean(tformTibFib.c);
        trab_tibia_fibulaF = trabF;
        
        %Apply same transformation to isloated fibula
        fibulaV = fibulaV * tformTibFib.T + mean(tformTibFib.c);
        
    else
        
        %Create the relevant data sets for the reference case
        %Tibia-fibula
        tibia_fibulaV = [tibiaV; fibulaV];
        tibia_fibulaF = [tibiaF(:,:,1); fibulaF(:,:,1) + length(tibiaV)];
        %Trabecular-tibia
        trab_tibiaV = trabV;
        trab_tibiaF = trabF;
        %Trabecular-tibia-fibula
        trab_tibia_fibulaV = trabV;
        trab_tibia_fibulaF = trabF;
        
    end
    
    %Store data in structure
    
    %Tibia models
    data.tibia.caseID{caseInd} = caseID{caseInd};
    data.tibia.F(:,:,caseInd) = tibiaF;
    data.tibia.V(:,:,caseInd) = tibiaV;
    
    %Fibula models
    data.fibula.caseID{caseInd} = caseID{caseInd};
    data.fibula.F(:,:,caseInd) = fibulaF;
    data.fibula.V(:,:,caseInd) = fibulaV;
    
    %Tibia-Fibula models
    data.tibia_fibula.caseID{caseInd} = caseID{caseInd};
    data.tibia_fibula.F(:,:,caseInd) = tibia_fibulaF;
    data.tibia_fibula.V(:,:,caseInd) = tibia_fibulaV;
    
    %Trabecular models
    %Two separate orientations given varying procrustes approaches
    %Trabecular-tibia
    data.trab_tibia.caseID{caseInd} = caseID{caseInd};
    data.trab_tibia.F(:,:,caseInd) = trab_tibiaF;
    data.trab_tibia.V(:,:,caseInd) = trab_tibiaV;
    %Trabecular-tibia-fibula
    data.trab_tibia_fibula.caseID{caseInd} = caseID{caseInd};
    data.trab_tibia_fibula.F(:,:,caseInd) = trab_tibia_fibulaF;
    data.trab_tibia_fibula.V(:,:,caseInd) = trab_tibia_fibulaV;
    
    %Export registered surfaces to STL
    
    %Create the stl structures
    %Tibia
    stlStructTib.solidNames = {['tibia-',caseNo]}; %names of parts
    stlStructTib.solidVertices = {tibiaV}; %Vertices
    stlStructTib.solidFaces = {tibiaF}; %Faces
    stlStructTib.solidNormals={[]};
    %Tibia-Fibula
    stlStructTibFib.solidNames = {['tibia-fibula-',caseNo]}; %names of parts
    stlStructTibFib.solidVertices = {tibia_fibulaV}; %Vertices
    stlStructTibFib.solidFaces = {tibia_fibulaF}; %Faces
    stlStructTibFib.solidNormals={[]};
    %Trabecular
    %Tibia version
    stlStructTrabTib.solidNames = {['trabecular-tibia',caseNo]}; %names of parts
    stlStructTrabTib.solidVertices = {trab_tibiaV}; %Vertices
    stlStructTrabTib.solidFaces = {trab_tibiaF}; %Faces
    stlStructTrabTib.solidNormals={[]};
    %Tibia-fibula version
    stlStructTrabTibFib.solidNames = {['trabecular-tibia-fibula',caseNo]}; %names of parts
    stlStructTrabTibFib.solidVertices = {trab_tibia_fibulaV}; %Vertices
    stlStructTrabTibFib.solidFaces = {trab_tibia_fibulaF}; %Faces
    stlStructTrabTibFib.solidNormals={[]};
    
    %Export STLs
    export_STL_txt(['..\..\ShapeModels\tibia\registered\',caseNo,'-tibia.stl'], stlStructTib);
    export_STL_txt(['..\..\ShapeModels\tibia-fibula\registered\',caseNo,'-tibia-fibula.stl'], stlStructTibFib);
    export_STL_txt(['..\..\ShapeModels\trabecular-tibia\registered\',caseNo,'-trabecular-tibia.stl'], stlStructTrabTib);
    export_STL_txt(['..\..\ShapeModels\trabecular-tibia-fibula\registered\',caseNo,'-trabecular-tibia-fibula.stl'], stlStructTrabTibFib);
    
    %Update waitbar
    waitbar(caseInd/length(caseID), wbar, ...
        ['Remeshed and aligned ',num2str(caseInd),' of ',num2str(length(caseID)),' surfaces...']);
    
    %Disp progress
    disp(['Case ',caseNo,' successfully remeshed and aligned.']);
    
    %Clear variables for re-loop
    clearvars -except visualiseRemesh caseID homeDir caseInd optionStruct_tib optionStruct_fib optCPD wbar data 
    
    %Navigate back to segmentation directory
    cd('..');
    
end

%Update and close waitbar
waitbar(1, wbar, 'Remeshed all surfaces. Finishing...');
close(wbar);

%Display finishing message
disp('FINISHED REGISTERING ALL SURFACES...');

%Save registered data points to .mat file
cd('..');
save('Data\registeredSurfaceDataPoints.mat', 'data');

%% Create shape models

%%%% TODO: include the below code in a conditional boolean as to whether to
%%%% register surfaces or load existing data --- registering takes a while
%%%% (~ 1-2 hours)

cd('..\data');
load('registeredSurfaceDataPoints.mat');

%% Tibia shape model

%Loop through case ID's and reshape data for PCA
for caseInd = 1:length(caseID)
    
    %Reshape data
    shapeModel.tibia.nodes(caseInd,:) = reshape(data.tibia.V(:,:,caseInd)', ...
        optionStruct_tib.nb_pts*3, [])';
    
end

%Run PCA to create shape model
[shapeModel.tibia.loadings, shapeModel.tibia.score, ...
    shapeModel.tibia.latent, ~, shapeModel.tibia.varExplained, ...
    shapeModel.tibia.mean] = pca(shapeModel.tibia.nodes);

%Create new faces variable to visualise
shapeModel.tibia.F = data.tibia.F(:,:,1);

%Reshape the mean points in the shape model to visualise
shapeModel.tibia.meanPoints = reshape(shapeModel.tibia.mean,...
    [3, length(shapeModel.tibia.mean)/3])';

% % % %Visualise the mean shape
% % % cFigure; hold on
% % % plotV(shapeModel.tibia.meanPoints,'b.');
% % % gpatch(shapeModel.tibia.F,shapeModel.tibia.meanPoints,'rw','k');
% % % axisGeom;

%Calculate the cumulative variance explained by the model
shapeModel.tibia.varExplainedSum = cumsum(shapeModel.tibia.varExplained);

%Identify PCs to retain that meet explanatory threshold
shapeModel.tibia.retainPCs = find(shapeModel.tibia.varExplainedSum > varExpThreshold,1);

%Reconstruct each case using the retained PCs
shapeModel.tibia.reconstructed = ...
    shapeModel.tibia.score(:,1:shapeModel.tibia.retainPCs) * ...
    shapeModel.tibia.loadings(:,1:shapeModel.tibia.retainPCs)';

% % % %Visualise a specific case
% % % %Specify case
% % % caseInd = input('Enter case index to visualise (1-30): ');
% % % %Reshape data structure and add the mean
% % % reconstructedV = reshape(shapeModel.tibia.reconstructed(caseInd,:) + shapeModel.tibia.mean, ...
% % %     [3, length(shapeModel.tibia.mean)/3])';
% % % %Visualise original and reconstructed surfaces
% % % cFigure; hold on;
% % % %Original
% % % gpatch(data.tibia.F(:,:,caseInd), data.tibia.V(:,:,caseInd),'gw','none',0.3);
% % % %Reconstructed
% % % gpatch(shapeModel.tibia.F, reconstructedV, 'rw','k');
% % % %Axes
% % % axisGeom;
% % % %Legend
% % % legend('Original', 'Reconstructed');

%Calculate the distance error between reconstructed and actual points across cases
for caseInd = 1:length(caseID)
    
    %Calculate point error distance
    shapeModel.tibia.pointErrorDist(caseInd,:) = ...
        distancePoints3d(data.tibia.V(:,:,caseInd), ...
        reshape((shapeModel.tibia.reconstructed(caseInd,:) + shapeModel.tibia.mean),...
        [3, length(shapeModel.tibia.mean)/3])');

    %Calculate the mean error
    shapeModel.tibia.pointErrorDistMean(caseInd,1) = ...
        mean(shapeModel.tibia.pointErrorDist(caseInd,:));
    
    %Calculate the peak error
    shapeModel.tibia.pointErrorDistMax(caseInd,1) = ...
        max(shapeModel.tibia.pointErrorDist(caseInd,:));
    
    %Convert distance error to colour scales for visualisation
    shapeModel.tibia.pointErrorDistColF(caseInd,:) = vertexToFaceMeasure(shapeModel.tibia.F, ...
        shapeModel.tibia.pointErrorDist(caseInd,:)');
    shapeModel.tibia.pointErrorDistColV(caseInd,:) = faceToVertexMeasure(shapeModel.tibia.F, ...
        data.tibia.V(:,:,caseInd), ...
        shapeModel.tibia.pointErrorDistColF(caseInd,:)');    
    
    %%%%% TODO: create a bone-red colour map for this??? Or a light to
    %%%%% black colormap???
    
    %%%%% TODO: does GIBBON offer that opportunity to smooth color,
    %%%%% extrapolate or something???
    
    %%%%% TODO: the view and rotating here could be looped with a changed
    %%%%% parameter to clean up code...
    
    %Create visualisation of error
    %Use subplots to create different perspectives
    cFigure; hold on;
    
    %Anterior view
    subplot(1,4,1);
    %Add surface
    hp = gpatch(shapeModel.tibia.F, ...
            reshape(shapeModel.tibia.reconstructed(caseInd,:) + shapeModel.tibia.mean, ...
            [3, length(shapeModel.tibia.mean)/3])', ...
            shapeModel.tibia.pointErrorDistColV(caseInd,:)', 'none', 1);
    %Interpolate colouring for smoothness
    hp.FaceColor = 'Interp'; colormap viridis
    %Set axis view
    axis equal; axis tight;
    view(0,90); rotate(hp,[0 1 0], -90);
    %Set axis parameters
    camlight headlight; axis off
    %Add title
    title('Anterior View', 'FontSize', 12);
    
    %Lateral view
    subplot(1,4,2);
    %Add surface
    hp = gpatch(shapeModel.tibia.F, ...
            reshape(shapeModel.tibia.reconstructed(caseInd,:) + shapeModel.tibia.mean, ...
            [3, length(shapeModel.tibia.mean)/3])', ...
            shapeModel.tibia.pointErrorDistColV(caseInd,:)', 'none', 1);
    %Interpolate colouring for smoothness
    hp.FaceColor = 'Interp'; colormap viridis
    %Set axis view
    axis equal; axis tight;
    view(0,90);
    %Set axis parameters
    camlight headlight; axis off
    %Add title
    title('Lateral View', 'FontSize', 12);
    
    %Posterior view
    subplot(1,4,3);
    %Add surface
    hp = gpatch(shapeModel.tibia.F, ...
            reshape(shapeModel.tibia.reconstructed(caseInd,:) + shapeModel.tibia.mean, ...
            [3, length(shapeModel.tibia.mean)/3])', ...
            shapeModel.tibia.pointErrorDistColV(caseInd,:)', 'none', 1);
    %Interpolate colouring for smoothness
    hp.FaceColor = 'Interp'; colormap viridis
    %Set axis view
    axis equal; axis tight;
    view(0,90); rotate(hp,[0 1 0], 90); 
    %Set axis parameters
    camlight headlight; axis off
    %Add title
    title('Posterior View', 'FontSize', 12);
    
    %Medial view
    subplot(1,4,4);
    %Add surface
    hp = gpatch(shapeModel.tibia.F, ...
            reshape(shapeModel.tibia.reconstructed(caseInd,:) + shapeModel.tibia.mean, ...
            [3, length(shapeModel.tibia.mean)/3])', ...
            shapeModel.tibia.pointErrorDistColV(caseInd,:)', 'none', 1);
    %Interpolate colouring for smoothness
    hp.FaceColor = 'Interp'; colormap viridis
    %Set axis view
    axis equal; axis tight;
    view(0,90); rotate(hp,[0 1 0], 180); 
    %Set axis parameters
    camlight headlight; axis off
    %Add colorbar on last view
    colorbar
    %Add title
    title('Medial View', 'FontSize', 12);
    
    %%%%% TODO: create a bone-red colour map for this??? Or a light to
    %%%%% black colormap???
        %%%% virids probably best option
    
    %%%%% TODO: the view and rotating here could be looped with a changed
    %%%%% parameter to clean up code...
    
    %%%%% TODO: finalise (e.g. fix colouring) and save plots...
    
end



%%%%% TODO: create figures to compare +/- PCs
%%%%% TODO: create animations to compare +/- PCs

%%%%% TODO: save shape model...

%%


%% Tibia-fibula shape model

%Loop through case ID's and reshape data for PCA
for caseInd = 1:length(caseID)
    
    %Reshape data
    shapeModel.tibia_fibula.nodes(caseInd,:) = reshape(data.tibia_fibula.V(:,:,caseInd)', ...
        (optionStruct_tib.nb_pts + optionStruct_fib.nb_pts)*3, [])';
    
end

%Run PCA to create shape model
[shapeModel.tibia_fibula.loadings, shapeModel.tibia_fibula.score, ...
    shapeModel.tibia_fibula.latent, ~, shapeModel.tibia_fibula.varExplained, ...
    shapeModel.tibia_fibula.mean] = pca(shapeModel.tibia_fibula.nodes);

%Create new faces variable to visualise
shapeModel.tibia_fibula.F = data.tibia_fibula.F(:,:,1);

%Reshape the mean points in the shape model to visualise
shapeModel.tibia_fibula.meanPoints = reshape(shapeModel.tibia_fibula.mean,...
    [3, length(shapeModel.tibia_fibula.mean)/3])';

% % % %Visualise the mean shape
% % % cFigure; hold on
% % % plotV(shapeModel.tibia_fibula.meanPoints,'b.');
% % % gpatch(shapeModel.tibia_fibula.F,shapeModel.tibia_fibula.meanPoints,'rw','k');
% % % axisGeom;

%Calculate the cumulative variance explained by the model
shapeModel.tibia_fibula.varExplainedSum = cumsum(shapeModel.tibia_fibula.varExplained);

%Identify PCs to retain that meet explanatory threshold
shapeModel.tibia_fibula.retainPCs = find(shapeModel.tibia_fibula.varExplainedSum > varExpThreshold,1);

%Reconstruct each case using the retained PCs
shapeModel.tibia_fibula.reconstructed = ...
    shapeModel.tibia_fibula.score(:,1:shapeModel.tibia_fibula.retainPCs) * ...
    shapeModel.tibia_fibula.loadings(:,1:shapeModel.tibia_fibula.retainPCs)';

% % % %Visualise a specific case
% % % %Specify case
% % % caseInd = input('Enter case index to visualise (1-30): ');
% % % %Reshape data structure and add the mean
% % % reconstructedV = reshape(shapeModel.tibia_fibula.reconstructed(caseInd,:) + shapeModel.tibia_fibula.mean, ...
% % %     [3, length(shapeModel.tibia_fibula.mean)/3])';
% % % %Visualise original and reconstructed surfaces
% % % cFigure; hold on;
% % % %Original
% % % gpatch(data.tibia_fibula.F(:,:,caseInd), data.tibia_fibula.V(:,:,caseInd),'gw','none',0.3);
% % % %Reconstructed
% % % gpatch(shapeModel.tibia_fibula.F, reconstructedV, 'rw','k');
% % % %Axes
% % % axisGeom;
% % % %Legend
% % % legend('Original', 'Reconstructed');

%Calculate the distance error between reconstructed and actual points across cases
for caseInd = 1:length(caseID)
    
    %Calculate point error distance
    shapeModel.tibia_fibula.pointErrorDist(caseInd,:) = ...
        distancePoints3d(data.tibia_fibula.V(:,:,caseInd), ...
        reshape((shapeModel.tibia_fibula.reconstructed(caseInd,:) + shapeModel.tibia_fibula.mean),...
        [3, length(shapeModel.tibia_fibula.mean)/3])');
    
    %Calculate the mean error
    shapeModel.tibia_fibula.pointErrorDistMean(caseInd,1) = ...
        mean(shapeModel.tibia_fibula.pointErrorDist(caseInd,:));
    
    %Calculate the peak error
    shapeModel.tibia_fibula.pointErrorDistMax(caseInd,1) = ...
        max(shapeModel.tibia_fibula.pointErrorDist(caseInd,:));
    
    %Convert distance error to colour faces for visualisation
    shapeModel.tibia_fibula.pointErrorDistCol(caseInd,:) = vertexToFaceMeasure(shapeModel.tibia_fibula.F, ...
        shapeModel.tibia_fibula.pointErrorDist(caseInd,:)');
    
    %%%%% TODO: add figure in...
    
    %%%%% TODO: finalise (e.g. fix colouring) and save plots...
    
end

%% Trabecular (tibia) shape model

%Loop through case ID's and reshape data for PCA
for caseInd = 1:length(caseID)
    
    %Reshape data
    shapeModel.trab_tibia.nodes(caseInd,:) = reshape(data.trab_tibia.V(:,:,caseInd)', ...
        optionStruct_tib.nb_pts*3, [])';
    
end

%Run PCA to create shape model
[shapeModel.trab_tibia.loadings, shapeModel.trab_tibia.score, ...
    shapeModel.trab_tibia.latent, ~, shapeModel.trab_tibia.varExplained, ...
    shapeModel.trab_tibia.mean] = pca(shapeModel.trab_tibia.nodes);

%Create new faces variable to visualise
shapeModel.trab_tibia.F = data.trab_tibia.F(:,:,1);

%Reshape the mean points in the shape model to visualise
shapeModel.trab_tibia.meanPoints = reshape(shapeModel.trab_tibia.mean,...
    [3, length(shapeModel.trab_tibia.mean)/3])';

% % % %Visualise the mean shape
% % % cFigure; hold on
% % % plotV(shapeModel.trab_tibia.meanPoints,'b.');
% % % gpatch(shapeModel.trab_tibia.F,shapeModel.trab_tibia.meanPoints,'rw','k');
% % % axisGeom;
% % % 
% % % %Visualise the mean shape of this relative to the tibia mean shape
% % % cFigure; hold on
% % % gpatch(shapeModel.tibia.F, shapeModel.tibia.meanPoints, 'gw', 'none', 0.3);
% % % gpatch(shapeModel.trab_tibia.F, shapeModel.trab_tibia.meanPoints, 'rw','k');
% % % axisGeom;

%Calculate the cumulative variance explained by the model
shapeModel.trab_tibia.varExplainedSum = cumsum(shapeModel.trab_tibia.varExplained);

%Identify PCs to retain that meet explanatory threshold
shapeModel.trab_tibia.retainPCs = find(shapeModel.trab_tibia.varExplainedSum > varExpThreshold,1);

%Reconstruct each case using the retained PCs
shapeModel.trab_tibia.reconstructed = ...
    shapeModel.trab_tibia.score(:,1:shapeModel.trab_tibia.retainPCs) * ...
    shapeModel.trab_tibia.loadings(:,1:shapeModel.trab_tibia.retainPCs)';

% % % %Visualise a specific case
% % % %Specify case
% % % caseInd = input('Enter case index to visualise (1-30): ');
% % % %Reshape data structure and add the mean
% % % reconstructedV = reshape(shapeModel.trab_tibia.reconstructed(caseInd,:) + shapeModel.trab_tibia.mean, ...
% % %     [3, length(shapeModel.trab_tibia.mean)/3])';
% % % %Visualise original and reconstructed surfaces
% % % cFigure; hold on;
% % % %Original
% % % gpatch(data.trab_tibia.F(:,:,caseInd), data.trab_tibia.V(:,:,caseInd),'gw','none',0.3);
% % % %Reconstructed
% % % gpatch(shapeModel.trab_tibia.F, reconstructedV, 'rw','k');
% % % %Axes
% % % axisGeom;
% % % %Legend
% % % legend('Original', 'Reconstructed');

%Calculate the distance error between reconstructed and actual points across cases
for caseInd = 1:length(caseID)
    
    %Calculate point error distance
    shapeModel.trab_tibia.pointErrorDist(caseInd,:) = ...
        distancePoints3d(data.trab_tibia.V(:,:,caseInd), ...
        reshape((shapeModel.trab_tibia.reconstructed(caseInd,:) + shapeModel.trab_tibia.mean),...
        [3, length(shapeModel.trab_tibia.mean)/3])');
    
    %Calculate the mean error
    shapeModel.trab_tibia.pointErrorDistMean(caseInd,1) = ...
        mean(shapeModel.trab_tibia.pointErrorDist(caseInd,:));
    
    %Calculate the peak error
    shapeModel.trab_tibia.pointErrorDistMax(caseInd,1) = ...
        max(shapeModel.trab_tibia.pointErrorDist(caseInd,:));
    
    %Convert distance error to colour faces for visualisation
    shapeModel.trab_tibia.pointErrorDistCol(caseInd,:) = vertexToFaceMeasure(shapeModel.trab_tibia.F, ...
        shapeModel.trab_tibia.pointErrorDist(caseInd,:)');
    
    %%%%% TODO: create a bone-red colour map for this??? Or a light to
    %%%%% black colormap???
    
    %%%%% TODO: does GIBBON offer that opportunity to smooth color,
    %%%%% extrapolate or something???
    
    %%%%% TODO: the view and rotating here could be looped with a changed
    %%%%% parameter to clean up code...
    
    %%%%% TODO: finalise (e.g. fix colouring) and save plots...
    
end

%%

%% Create animation of shape change for PC

%%%% TODO: use shift to vertex colouring and interpolate in animations...

%%%% TODO: embed into function...

%Set PC
PC = 4;

%Set number of PCs
reconPCs = length(shapeModel.trab_tibia.varExplained);

%Reshape mean
meanPts = reshape(shapeModel.trab_tibia.mean, [3, length(shapeModel.trab_tibia.mean)/3])';

%Calculate mean and standard deviation values for scores
% % % meanScore = mean(shapeModel.score); %basically zero
sdScore = std(shapeModel.trab_tibia.score);

%Create a variable to work through 0.1 increments from -3 to +3 SD
calcSD = linspace(0.1, 5, diff([0,5])*10);

%Calculate the points at each increment
for nSD = 1:length(calcSD)
    
    %Reconstruct the points at the current SD interval

    %Set the simulated scores array with the standard deviation added/subtracted
    %Plus SD
    simScorePlus = zeros(1,reconPCs);
    simScorePlus(PC) = simScorePlus(PC) + (sdScore(PC) * calcSD(nSD));
    %Minus SD
    simScoreMinus = zeros(1,reconPCs);
    simScoreMinus(PC) = simScoreMinus(PC) + (sdScore(PC) * -calcSD(nSD));
    
    %Reconstruct the points, add the mean and reshape to 3D
    plusPts = reshape((simScorePlus * shapeModel.trab_tibia.loadings(:,1:reconPCs)') + shapeModel.trab_tibia.mean, ...
        [3, length(shapeModel.trab_tibia.mean)/3])';
    minusPts = reshape((simScoreMinus * shapeModel.trab_tibia.loadings(:,1:reconPCs)') + shapeModel.trab_tibia.mean, ...
        [3, length(shapeModel.trab_tibia.mean)/3])';
        
    %Set into structure
    posV(:,:,nSD) = plusPts;
    negV(:,:,nSD) = minusPts;
    
    %Calculate distances between points for heat mapping
    pcPlusDist = distancePoints3d(meanPts, plusPts);
    pcMinusDist = distancePoints3d(meanPts, minusPts);

    %Convert distances to face measures
%     pcPlusC = vertexToFaceMeasure(data.trab_tibia.F(:,:,1), pcPlusDist);
%     pcMinusC = vertexToFaceMeasure(data.trab_tibia.F(:,:,1), pcMinusDist);
    pcPlusC = vertexToFaceMeasure(shapeModel.trab_tibia.F, pcPlusDist);
    pcMinusC = vertexToFaceMeasure(shapeModel.trab_tibia.F, pcMinusDist);

    %Set into structure
    posC(:,nSD) = pcPlusC;
    negC(:,nSD) = pcMinusC; 
    
end

%Create the basic view to store graphics and initiate animation
hf = cFigure;
%Positive SD
subplot(1,2,1);
title(['Plus SD for PC',num2str(PC)]);
% hp_pos = gpatch(data.trab_tibia.F(:,:,1), meanPts, posC(:,1), 'none', 1);
hp_pos = gpatch(shapeModel.trab_tibia.F, meanPts, posC(:,1), 'none', 1);
axisGeom; colormap viridis; camlight headlight
axis(axisLim([posV;negV])); %set axis limits based on values
axis off
%Negative SD
subplot(1,2,2);
title(['Minus SD for PC',num2str(PC)]);
% hp_neg = gpatch(data.trab_tibia.F(:,:,1), meanPts, negC(:,1), 'none', 1);
hp_neg = gpatch(shapeModel.trab_tibia.F(:,:,1), meanPts, negC(:,1), 'none', 1);
axisGeom; colormap viridis; camlight headlight
axis(axisLim([posV;negV])); %set axis limits based on values
axis off

%Set up animation
animStruct.Time = calcSD;
%Loop through SD points
for nSD = 1:length(calcSD)
    
    %Set entries into animation structure
    animStruct.Handles{nSD} = [hp_pos hp_pos hp_neg hp_neg]; %Handles of objects to animate
    animStruct.Props{nSD} = {'Vertices', 'CData', 'Vertices', 'CData'}; %Properties of objects to animate
    animStruct.Set{nSD} = {posV(:,:,nSD), posC(:,nSD), negV(:,:,nSD), negC(:,nSD)}; %Property values for to set in order to animate
    
end

%Animate figure
anim8(hf, animStruct);

% % % %Create a static version
% % % cFigure; hold on;
% % % %Mean shape - grey
% % % hp_m = gpatch(shapeModel.trab_tibia.F, meanPts, [200/255 200/255 200/255], 'none', 0.7);
% % % %Plus SD
% % % hp_plusSD = gpatch(shapeModel.trab_tibia.F, posV(:,:,length(calcSD)), 'gw', 'none', 0.7);
% % % %Minus SD
% % % hp_minusSD = gpatch(shapeModel.trab_tibia.F, negV(:,:,length(calcSD)), 'rw', 'none', 0.7);
% % % %Axes
% % % axisGeom; camlight headlight

% % % exportGifAnim8();

%%%% TODO: could provide alternative views across subplots in anim8 figure?

%%%% NOTE: for PC1 of tib-fib there's a bit of noise at larger SDs as the
%%%% nodes get squished together or pulled apart --- unsure of how to
%%%% address that??? It occurs across other components too --- it's like a
%%%% distortion --- but that might be as I'm testing it with only 10
%%%% surfaces...
%%%%
%%%% Outside of this, the tibia-fibula model seems to be picking up typical
%%%% shape variation versus just random noise...
%%%%
%%%% Doesn't get much better necessarily with more cases in the shape model
%%%% - but it's perhaps a little more noticeable with the color mapping, as
%%%% coincident points have more color variation


%% Test out regression model between tibia and trabecular shapes

%Regression model for trabecular PC

%%%% TODO: cross validation k-fold approach (70:30 split) to train
%%%% regression model???

%Create prediction matrix
%Use retained PCs as test
X = shapeModel.tibia.score(:,1:shapeModel.tibia.retainPCs);
% % % X = zscore(shapeModel.tibia.score(:,1:shapeModel.tibia.retainPCs));

%Loop through trabecular PCs to retain
for predictPC = 1:shapeModel.trab_tibia.retainPCs

    %Create output matrix
    y = shapeModel.trab_tibia.score(:,predictPC);
    % % % y = zscore(shapeModel.trab_tibia.score(:,predictPC));

    %Would removing outliers > 1 or 2 SD assist in creating a better model???
    %Going to < 1 SD removes too many data points
    %For predicting PC2 for the trabecular it does improve the predictive model
    %a bit (R2 from 0.4 to 0.6). Going down to < 1 makes it worse again
    %Same approach actually makes PC3 slightly worse
    %Makes PC4 a bit worse
    % X = X(find(abs(y) < 2),:);
    % y = y(find(abs(y) < 2),:);

    %Create linear model
    regMdl.trab_tibia{predictPC} = fitlm(X,y);
    %%%% Gives pretty good R-squared value for PC1 (R-squared > 0.9)...
    %%%% R-squared for PC2 not as good (R-squared ~0.4)...outlier...
    %%%% PC3 a bit better (R-squared ~0.63)...less obvious outliers...
    %%%% PC4 back down to R-squared of ~0.32...
    %%%%
    %%%% Will be interesting to see how this comes out in overall shape error,
    %%%% the ability to predict the highest proportion of variance component
    %%%% quite well could assist in minimising error...
    
    %%%% Cross-validation option
% % %     %Run model
% % %     regMdl_CV = fitrlinear(X,y,'CrossVal','on','KFold',10, ...
% % %         'ObservationsIn', 'rows');
% % %     %Predict responses for out of fold observations
% % %     yHat_OoF = kfoldPredict(regMdl_CV);
    

    %%%% Manual leave-one-out cross-validation approach
    %%%%%% Essentially you could do the same linear regression approach
    %%%%%% currently being used --- but randomly leave a sub-section of
    %%%%%% data out. You would then predict those that weren't in there and
    %%%%%% each time you'd get an error
    %%%%%% e.g.
% % %     regMdl_leaveOut = fitlm(X(1:27,:),y(1:27,:));
% % %     yPred_leaveOut = predict(regMdl_leaveOut, X(28:30,:));
% % %     leaveOut_mse = mean(sqrt((yPred_leaveOut - y(28:30)).^2));
    %%%%%% The mean squared error is somewhat pointless for an individual
    %%%%%% component though --- as the units are dimensionless --- probably
    %%%%%% better to run all PC models using the same leave-out criteria,
    %%%%%% reconstruct the surfaces and calculate the mean and max point
    %%%%%% error estimates...

    %Output predicted values
    yPred(:,predictPC) = predict(regMdl.trab_tibia{predictPC}, X);

% % %     %Calculate the root mean square error between the two
% % %     rmse = sqrt(sum((y - yPred).^2) / length(y));
    
end

%Reconstruct each case using the predicted PCs
reconstructPred = ...
    yPred(:,1:shapeModel.trab_tibia.retainPCs) * ...
    shapeModel.trab_tibia.loadings(:,1:shapeModel.trab_tibia.retainPCs)';

%Visualise a specific case
%Specify case
caseInd = input('Enter case index to visualise (1-30): ');
%Reshape data structure and add the mean
%Original reconstruction from shape model
reconstructedV = reshape(shapeModel.trab_tibia.reconstructed(caseInd,:) + shapeModel.trab_tibia.mean, ...
    [3, length(shapeModel.trab_tibia.mean)/3])';
%Predicted reconstruction from regression model
predictedV = reshape(reconstructPred(caseInd,:) + shapeModel.trab_tibia.mean, ...
    [3, length(shapeModel.trab_tibia.mean)/3])';

%Visualise original, reconstructed and predicted surfaces
cFigure; hold on;

%Original vs. reconstructed from shape model
subplot(1,2,1); hold on;
gpatch(data.trab_tibia.F(:,:,caseInd), data.trab_tibia.V(:,:,caseInd),'gw','none',0.3);
gpatch(shapeModel.trab_tibia.F, reconstructedV, 'rw','k');
%Axes
axisGeom;
%Title
title('Original vs. Shape Model Reconstruction', 'FontSize', 14);

%Original vs. predicted from regression model
subplot(1,2,2); hold on;
gpatch(data.trab_tibia.F(:,:,caseInd), data.trab_tibia.V(:,:,caseInd),'gw','none',0.3);
gpatch(shapeModel.trab_tibia.F, predictedV, 'rw','k');
%Axes
axisGeom;
%Title
title('Original vs. Regression Model Reconstruction', 'FontSize', 14);

%%%% TODO: compare errors to both original surface and reconstructed
%%%% surface from the same retained components




