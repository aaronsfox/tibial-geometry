%% This script provides the steps to register the 3D bone surfaces in 
%  preparation for creating the statistical shape model. Components within
%  this function are time consuming, and hence options are provided to load
%  the pre-existing data vs. re-running the processes.
%
%  See the README.MD in the base level folder for more descriptive details
%  on this process and the source data.
%
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

%% Inputs

%Flag to run the registration process
%Defaulted to false as this is a time consuming process (~1-2 hours)
registerSurfaces = false;

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
optionStruct_tib.nb_pts = 3500; %Set desired number of points
optionStruct_tib.disp_on = 0; % Turn off command window text display
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

%Set threshold for keeping PCs of models (i.e. % variance explained)
varExpThreshold = 95;

%% Remesh and align surfaces

%All surface meshes need to be aligned to a single target surface. For this
%we use a case that is closest to the mean surface in line with existing
% work (i.e. Bruce et al. 2021, Computer Methods Biomech Biomed Eng, doi:
%10.1080/10255842.2021.1985111). Some testing revealed that the best use
%case for this was case-147211, hence we ensure this is processed first and
%all subsequent meshes are aligned to this

%Check flag
if registerSurfaces

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
        tibPoints = [{'LC'},{'LM'},{'MC'},{'MM'}];
        for pp = 1:length(tibPoints)
            tree = xml_read([tibPoints{pp},'.txt']);
            landmarks.(char(tree.Point.Name)) = tree.Point.Coordinate;
            clear tree 
        end
        clear pp

        %Create new landmarks
        landmarks.IM = midPoint3d(landmarks.LM,landmarks.MM);
        landmarks.IC = midPoint3d(landmarks.LC,landmarks.MC);

        %Create and align tibial coordinate system
        [tibiaF, tibiaV, fibulaF, fibulaV, trabF, trabV] = ...
            rotateSurfacesGlobal(landmarks, tibiaF, tibiaV, fibulaF, fibulaV, trabF, trabV);

        %If not the target mesh, align the surfaces
        if caseInd > 1

            %Tibia

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp(['APPLYING CPD ALGORITHM TO TIBIA OF CASE-',caseNo,'...']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %Perform non-rigid registration using CPD to align surfaces
            [cpdTformTib] = cpd_register(data.tibia.V(:,:,1), tibiaV, optCPD);

            %Identify the matching points in the target mesh against the
            %registration to identify the corresponding point indices
            regSortIdxTib = knnsearch(cpdTformTib.Y, data.tibia.V(:,:,1));

            %Sort the registered points so they align with the reference mesh
            tibiaV_sorted = tibiaV(regSortIdxTib,:);

            %Visualise to confirm that points now correspond. This can be done
            %by using the target mesh faces against the registered moving
            %points. The surface meshes should essentially come out looking the
            %same (i.e. overlapped) despite using different faces
            cFigure; hold on;
            subplot(1,3,1); hold on;
            gpatch(tibiaF, tibiaV, 'rw', 'none', 0.3);
            gpatch(data.tibia.F(:,:,1), tibiaV_sorted, 'gw', 'k');
            title('Registered Tibia');
            axisGeom;;

            %Allocate to original variables
            tibiaF = data.tibia.F(:,:,1); %use target mesh faces
            tibiaV = tibiaV_sorted; %use sorted registered points

            %Trabecular

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp(['APPLYING CPD ALGORITHM TO TRABECULAR OF CASE-',caseNo,'...']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %Perform non-rigid registration using CPD to align surfaces
            [cpdTformTrab] = cpd_register(data.trab_tibia.V(:,:,1), trabV, optCPD);

            %Identify the matching points in the target mesh against the
            %registration to identify the corresponding point indices
            regSortIdxTrab = knnsearch(cpdTformTrab.Y, data.trab_tibia.V(:,:,1));

            %Sort the registered points so they align with the shape model mean points
            trabV_sorted = trabV(regSortIdxTrab,:);

            %Visualise to confirm that points now correspond. This can be done
            %by using the target mesh faces against the registered moving
            %points. The surface meshes should essentially come out looking the
            %same (i.e. overlapped) despite using different faces
            subplot(1,3,2); hold on;
            gpatch(trabF, trabV, 'rw', 'none', 0.3);
            gpatch(data.trab_tibia.F(:,:,1), trabV_sorted, 'gw', 'k');
            title('Registered Trabecular');
            axisGeom;
            
            %Allocate to original variables
            trabF = data.trab_tibia.F(:,:,1); %use target mesh faces
            trabV = trabV_sorted; %use sorted registered points

            %Fibula

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp(['APPLYING CPD ALGORITHM TO FIBULA OF CASE-',caseNo,'...']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %Perform non-rigid registration using CPD to align surfaces
            [cpdTformFib] = cpd_register(data.fibula.V(:,:,1), fibulaV, optCPD);

            %Identify the matching points in the target mesh against the
            %registration to identify the corresponding point indices
            regSortIdxFib = knnsearch(cpdTformFib.Y, data.fibula.V(:,:,1));

            %Sort the registered points so they align with the shape model mean points
            fibulaV_sorted = fibulaV(regSortIdxFib,:);

            %Visualise to confirm that points now correspond. This can be done
            %by using the target mesh faces against the registered moving
            %points. The surface meshes should essentially come out looking the
            %same (i.e. overlapped) despite using different faces
            subplot(1,3,3); hold on;
            gpatch(fibulaF, fibulaV, 'rw', 'none', 0.3);
            gpatch(data.fibula.F(:,:,1), fibulaV_sorted, 'gw', 'k');
            title('Registered Fibula');
            axisGeom;

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
        clearvars -except caseID registerSurfaces homeDir caseInd optionStruct_tib optionStruct_fib optCPD wbar data 

        %Navigate back to segmentation directory
        cd('..');

    end

    %Update and close waitbar
    waitbar(1, wbar, 'Remeshed all surfaces. Finishing...');
    close(wbar);

    %Display finishing message
    disp('FINISHED REGISTERING ALL SURFACES...');

    %Save registered data points to .mat file
    cd('..\Data');
    save('registeredSurfaceDataPoints.mat', 'data');
    
else
    
    %Load the existing registered surface dataset
    cd('..\Data');
    load('registeredSurfaceDataPoints.mat');
    
end

%% Create shape models

%Navigate to shape models directory
cd('..\ShapeModels');

%Set array to rotate views of surfaces later in reconstruction process
surfaceRot = [-90, 0, 90, 180];

%Set list to label reconstruction subplots with
viewLabel = [{'Anterior'}, {'Lateral'}, {'Posterior'}, {'Medial'}];

%% Tibia shape model

%Navigate into shape model directory
cd('tibia');

%Loop through case ID's and reshape data for PCA
for caseInd = 1:length(caseID)
    
    %Reshape data
    tibiaShapeModel.nodes(caseInd,:) = reshape(data.tibia.V(:,:,caseInd)', ...
        optionStruct_tib.nb_pts*3, [])';
    
end

%Run PCA to create shape model
[tibiaShapeModel.loadings, tibiaShapeModel.score, ...
    tibiaShapeModel.latent, ~, tibiaShapeModel.varExplained, ...
    tibiaShapeModel.mean] = pca(tibiaShapeModel.nodes);

%Create new faces variable to visualise
tibiaShapeModel.F = data.tibia.F(:,:,1);

%Reshape the mean points in the shape model to visualise
tibiaShapeModel.meanPoints = reshape(tibiaShapeModel.mean,...
    [3, length(tibiaShapeModel.mean)/3])';

%Calculate the cumulative variance explained by the model
tibiaShapeModel.varExplainedSum = cumsum(tibiaShapeModel.varExplained);

%Identify PCs to retain that meet explanatory threshold
tibiaShapeModel.retainPCs = find(tibiaShapeModel.varExplainedSum > varExpThreshold,1);

%Reconstruct each case using the retained PCs
tibiaShapeModel.reconstructed = ...
    tibiaShapeModel.score(:,1:tibiaShapeModel.retainPCs) * ...
    tibiaShapeModel.loadings(:,1:tibiaShapeModel.retainPCs)';

%Calculate the distance error between reconstructed and actual points across cases
for caseInd = 1:length(caseID)
    
    %Calculate point error distance
    tibiaShapeModel.pointErrorDist(caseInd,:) = ...
        distancePoints3d(data.tibia.V(:,:,caseInd), ...
        reshape((tibiaShapeModel.reconstructed(caseInd,:) + tibiaShapeModel.mean),...
        [3, length(tibiaShapeModel.mean)/3])');

    %Calculate the mean error
    tibiaShapeModel.pointErrorDistMean(caseInd,1) = ...
        mean(tibiaShapeModel.pointErrorDist(caseInd,:));
    
    %Calculate the peak error
    tibiaShapeModel.pointErrorDistMax(caseInd,1) = ...
        max(tibiaShapeModel.pointErrorDist(caseInd,:));
    
    %Convert distance error to colour scales for visualisation
    tibiaShapeModel.pointErrorDistColF(caseInd,:) = vertexToFaceMeasure(tibiaShapeModel.F, ...
        tibiaShapeModel.pointErrorDist(caseInd,:)');
    tibiaShapeModel.pointErrorDistColV(caseInd,:) = faceToVertexMeasure(tibiaShapeModel.F, ...
        data.tibia.V(:,:,caseInd), ...
        tibiaShapeModel.pointErrorDistColF(caseInd,:)');    
    
    %Create visualisation of error
    %Use subplots to create different perspectives
    cFigure; hold on;    
    %Loop through four views to create subplot
    for viewNo = 1:4    
        %Create subplot for current view
        subplot(1,4,viewNo);
        %Add surface
        hp = gpatch(tibiaShapeModel.F, ...
                reshape(tibiaShapeModel.reconstructed(caseInd,:) + tibiaShapeModel.mean, ...
                [3, length(tibiaShapeModel.mean)/3])', ...
                tibiaShapeModel.pointErrorDistColV(caseInd,:)', 'none', 1);
        %Interpolate colouring for smoothness
        hp.FaceColor = 'Interp'; colormap viridis
        %Set axis view
        axis equal; axis tight;
        view(0,90); rotate(hp,[0 1 0], surfaceRot(viewNo));
        %Set axis parameters
        camlight headlight; axis off
        %Add colorbar on last view
        if viewNo == 4
            colorbar
        end
        %Add title
        title([viewLabel{viewNo},' View'], 'FontSize', 12);
    end
    
    %Export figure
    export_fig(['figures\reconstructionErrors\',caseID{caseInd},'_reconstructionErrorMap.png'],'-m1');
	close
    
end

%This function can be used to examine the shape effect of a principal
%component. It creates an interactive animation of increasing/decreasing
%the score of a specific principal component by a factor of the standard
%deviation (i.e. the fairly standardised method for interpreting principal
%components). The below code is defaulting to visualising PC1 over a -3 to
%+3 standard deviation range.
% % % animatePrincipalComponent(tibiaShapeModel, 1, 3);

%After running the above animate function there is an option to export the
%gif using further GIBBON functionality using the command below. Refer to
%GIBBON documentation for further details on anim8 functionality.
% % % exportGifAnim8();

%Save the shape model data
save('tibiaShapeModel.mat', 'tibiaShapeModel');

%Return to base shape model directory
cd('..');

%% Tibia-fibula shape model

%Navigate into shape model directory
cd('tibia-fibula');

%Loop through case ID's and reshape data for PCA
for caseInd = 1:length(caseID)
    
    %Reshape data
    tibiaFibulaShapeModel.nodes(caseInd,:) = reshape(data.tibia_fibula.V(:,:,caseInd)', ...
        (optionStruct_tib.nb_pts + optionStruct_fib.nb_pts)*3, [])';
    
end

%Run PCA to create shape model
[tibiaFibulaShapeModel.loadings, tibiaFibulaShapeModel.score, ...
    tibiaFibulaShapeModel.latent, ~, tibiaFibulaShapeModel.varExplained, ...
    tibiaFibulaShapeModel.mean] = pca(tibiaFibulaShapeModel.nodes);

%Create new faces variable to visualise
tibiaFibulaShapeModel.F = data.tibia_fibula.F(:,:,1);

%Reshape the mean points in the shape model to visualise
tibiaFibulaShapeModel.meanPoints = reshape(tibiaFibulaShapeModel.mean,...
    [3, length(tibiaFibulaShapeModel.mean)/3])';

%Calculate the cumulative variance explained by the model
tibiaFibulaShapeModel.varExplainedSum = cumsum(tibiaFibulaShapeModel.varExplained);

%Identify PCs to retain that meet explanatory threshold
tibiaFibulaShapeModel.retainPCs = find(tibiaFibulaShapeModel.varExplainedSum > varExpThreshold,1);

%Reconstruct each case using the retained PCs
tibiaFibulaShapeModel.reconstructed = ...
    tibiaFibulaShapeModel.score(:,1:tibiaFibulaShapeModel.retainPCs) * ...
    tibiaFibulaShapeModel.loadings(:,1:tibiaFibulaShapeModel.retainPCs)';

%Calculate the distance error between reconstructed and actual points across cases
for caseInd = 1:length(caseID)
    
    %Calculate point error distance
    tibiaFibulaShapeModel.pointErrorDist(caseInd,:) = ...
        distancePoints3d(data.tibia_fibula.V(:,:,caseInd), ...
        reshape((tibiaFibulaShapeModel.reconstructed(caseInd,:) + tibiaFibulaShapeModel.mean),...
        [3, length(tibiaFibulaShapeModel.mean)/3])');

    %Calculate the mean error
    tibiaFibulaShapeModel.pointErrorDistMean(caseInd,1) = ...
        mean(tibiaFibulaShapeModel.pointErrorDist(caseInd,:));
    
    %Calculate the peak error
    tibiaFibulaShapeModel.pointErrorDistMax(caseInd,1) = ...
        max(tibiaFibulaShapeModel.pointErrorDist(caseInd,:));
    
    %Convert distance error to colour scales for visualisation
    tibiaFibulaShapeModel.pointErrorDistColF(caseInd,:) = vertexToFaceMeasure(tibiaFibulaShapeModel.F, ...
        tibiaFibulaShapeModel.pointErrorDist(caseInd,:)');
    tibiaFibulaShapeModel.pointErrorDistColV(caseInd,:) = faceToVertexMeasure(tibiaFibulaShapeModel.F, ...
        data.tibia_fibula.V(:,:,caseInd), ...
        tibiaFibulaShapeModel.pointErrorDistColF(caseInd,:)'); 
    
    %Create visualisation of error
    %Use subplots to create different perspectives
    cFigure; hold on;    
    %Loop through four views to create subplot
    for viewNo = 1:4    
        %Create subplot for current view
        subplot(1,4,viewNo);
        %Add surface
        hp = gpatch(tibiaFibulaShapeModel.F, ...
                reshape(tibiaFibulaShapeModel.reconstructed(caseInd,:) + tibiaFibulaShapeModel.mean, ...
                [3, length(tibiaFibulaShapeModel.mean)/3])', ...
                tibiaFibulaShapeModel.pointErrorDistColV(caseInd,:)', 'none', 1);
        %Interpolate colouring for smoothness
        hp.FaceColor = 'Interp'; colormap viridis
        %Set axis view
        axis equal; axis tight;
        view(0,90); rotate(hp,[0 1 0], surfaceRot(viewNo));
        %Set axis parameters
        camlight headlight; axis off
        %Add colorbar on last view
        if viewNo == 4
            colorbar
        end
        %Add title
        title([viewLabel{viewNo},' View'], 'FontSize', 12);
    end
    
    %Export figure
    export_fig(['figures\reconstructionErrors\',caseID{caseInd},'_reconstructionErrorMap.png'],'-m1');
	close
    
end

%This function can be used to examine the shape effect of a principal
%component. It creates an interactive animation of increasing/decreasing
%the score of a specific principal component by a factor of the standard
%deviation (i.e. the fairly standardised method for interpreting principal
%components). The below code is defaulting to visualising PC1 over a -3 to
%+3 standard deviation range.
% % % animatePrincipalComponent(tibiaFibulaShapeModel, 1, 3);

%After running the above animate function there is an option to export the
%gif using further GIBBON functionality using the command below. Refer to
%GIBBON documentation for further details on anim8 functionality.
% % % exportGifAnim8();

%Save the shape model data
save('tibiaFibulaShapeModel.mat', 'tibiaFibulaShapeModel');

%Return to base shape model directory
cd('..');

%% Trabecular (tibia) shape model

%Navigate into shape model directory
cd('trabecular-tibia');

%Loop through case ID's and reshape data for PCA
for caseInd = 1:length(caseID)
    
    %Reshape data
    trabShapeModel.nodes(caseInd,:) = reshape(data.trab_tibia.V(:,:,caseInd)', ...
        optionStruct_tib.nb_pts*3, [])';
    
end

%Run PCA to create shape model
[trabShapeModel.loadings, trabShapeModel.score, ...
    trabShapeModel.latent, ~, trabShapeModel.varExplained, ...
    trabShapeModel.mean] = pca(trabShapeModel.nodes);

%Create new faces variable to visualise
trabShapeModel.F = data.trab_tibia.F(:,:,1);

%Reshape the mean points in the shape model to visualise
trabShapeModel.meanPoints = reshape(trabShapeModel.mean,...
    [3, length(trabShapeModel.mean)/3])';

%Calculate the cumulative variance explained by the model
trabShapeModel.varExplainedSum = cumsum(trabShapeModel.varExplained);

%Identify PCs to retain that meet explanatory threshold
trabShapeModel.retainPCs = find(trabShapeModel.varExplainedSum > varExpThreshold,1);

%Reconstruct each case using the retained PCs
trabShapeModel.reconstructed = ...
    trabShapeModel.score(:,1:trabShapeModel.retainPCs) * ...
    trabShapeModel.loadings(:,1:trabShapeModel.retainPCs)';

%Calculate the distance error between reconstructed and actual points across cases
for caseInd = 1:length(caseID)
    
    %Calculate point error distance
    trabShapeModel.pointErrorDist(caseInd,:) = ...
        distancePoints3d(data.trab_tibia.V(:,:,caseInd), ...
        reshape((trabShapeModel.reconstructed(caseInd,:) + trabShapeModel.mean),...
        [3, length(trabShapeModel.mean)/3])');

    %Calculate the mean error
    trabShapeModel.pointErrorDistMean(caseInd,1) = ...
        mean(trabShapeModel.pointErrorDist(caseInd,:));
    
    %Calculate the peak error
    trabShapeModel.pointErrorDistMax(caseInd,1) = ...
        max(trabShapeModel.pointErrorDist(caseInd,:));
    
    %Convert distance error to colour scales for visualisation
    trabShapeModel.pointErrorDistColF(caseInd,:) = vertexToFaceMeasure(trabShapeModel.F, ...
        trabShapeModel.pointErrorDist(caseInd,:)');
    trabShapeModel.pointErrorDistColV(caseInd,:) = faceToVertexMeasure(trabShapeModel.F, ...
        data.trab_tibia.V(:,:,caseInd), ...
        trabShapeModel.pointErrorDistColF(caseInd,:)'); 
    
    %Create visualisation of error
    %Use subplots to create different perspectives
    cFigure; hold on;    
    %Loop through four views to create subplot
    for viewNo = 1:4    
        %Create subplot for current view
        subplot(1,4,viewNo);
        %Add surface
        hp = gpatch(trabShapeModel.F, ...
                reshape(trabShapeModel.reconstructed(caseInd,:) + trabShapeModel.mean, ...
                [3, length(trabShapeModel.mean)/3])', ...
                trabShapeModel.pointErrorDistColV(caseInd,:)', 'none', 1);
        %Interpolate colouring for smoothness
        hp.FaceColor = 'Interp'; colormap viridis
        %Set axis view
        axis equal; axis tight;
        view(0,90); rotate(hp,[0 1 0], surfaceRot(viewNo));
        %Set axis parameters
        camlight headlight; axis off
        %Add colorbar on last view
        if viewNo == 4
            colorbar
        end
        %Add title
        title([viewLabel{viewNo},' View'], 'FontSize', 12);
    end
    
    %Export figure
    export_fig(['figures\reconstructionErrors\',caseID{caseInd},'_reconstructionErrorMap.png'],'-m1');
	close
    
end

%This function can be used to examine the shape effect of a principal
%component. It creates an interactive animation of increasing/decreasing
%the score of a specific principal component by a factor of the standard
%deviation (i.e. the fairly standardised method for interpreting principal
%components). The below code is defaulting to visualising PC1 over a -3 to
%+3 standard deviation range.
% % % animatePrincipalComponent(trabShapeModel, 1, 3);

%After running the above animate function there is an option to export the
%gif using further GIBBON functionality using the command below. Refer to
%GIBBON documentation for further details on anim8 functionality.
% % % exportGifAnim8();

%Save the shape model data
save('trabShapeModel.mat', 'trabShapeModel');

%Return to base shape model directory
cd('..');

%% ----- end of createShapeModel.m ----- %%