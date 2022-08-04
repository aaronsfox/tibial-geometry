%% This script provides a demonstration of how to generate the tibia-fibula
%  surface from a simple set of landmarks and the statistical shape model,
%  and provides an evaluation of the error associated with this across a
%  set of cases.
%
%  See the README.MD in this folder for more descriptive details on this
%  process.
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

%% Options

%Flag to run all reconstructions and calculations or load in results
%Default to false as this takes time
runReconstructions = false;

%%%%%%% TODO: include this conditional throughout code

%% Set-up

%Set home directory
homeDir = pwd;

%Add supplementary code paths
addpath(genpath('supplementary'));

%Set options for remeshing
optionStruct_tib.nb_pts = 3500; %Set desired number of points
optionStruct_tib.disp_on = 0; % Turn off command window text display
optionStruct_fib.nb_pts = 2000; %Set desired number of points
optionStruct_fib.disp_on = 0; % Turn off command window text display

%Set options for CPD algorithm - non rigid
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
optCPD.corresp = 1;         % estimate correspondence

%Set options for CPD algorithm - rigid
optCPD_rig.method = 'rigid';    % use rigid registration
optCPD_rig.viz = 0;             % don't visualise
optCPD_rig.outliers = 0;        % don't account for outliers
optCPD_rig.fgt = 0;             % do not use FGT (default)
optCPD_rig.normalize = 1;       % normalize to unit variance and zero mean before registering (default)
optCPD_rig.corresp = 0;         % compute correspondence vector at the end of registration (not being estimated by default)
optCPD_rig.max_it = 100;        % max number of iterations
optCPD_rig.tol = 1e-4;          % tolerance
optCPD_rig.corresp = 0;         % estimate correspondence
optCPD_rig.scale = 0;           % turn off scaling

%Set array to rotate views of surfaces later in reconstruction process
surfaceRot = [-90, 0, 90, 180];

%Set list to label reconstruction subplots with
viewLabel = [{'Anterior'}, {'Lateral'}, {'Posterior'}, {'Medial'}];

%Turn warnings off
warning off

%Set list of palpable landmarks
palpableLandmarks =[{'MC'}; {'TT'}; {'T25'}; {'T50'}; {'T75'};
    {'MM'}; {'FH'}; {'F75'}; {'LM'}];

%% Load in and process mean data of shape model

%Load the tibia-fibula shape model required
load('..\..\ShapeModels\tibia-fibula\tibiaFibulaShapeModel.mat');

%Load shape model landmarks
shapeModelLandmarks = readtable('..\..\ShapeModels\tibia-fibula\tibiaFibulaShapeModel_mean.csv');
shapeModelLandmarks.Properties.VariableNames = {'landmark' 'X' 'Y' 'Z'};

%Extract landmarks into structure as variables
for landmarkNo = 1:length(shapeModelLandmarks.landmark)
    modelLandmarks.(shapeModelLandmarks.landmark{landmarkNo}) = ...
        [shapeModelLandmarks.X(landmarkNo), ...
        shapeModelLandmarks.Y(landmarkNo), ...
        shapeModelLandmarks.Z(landmarkNo)];    
end

%Create new landmarks for mean model
modelLandmarks.IM = midPoint3d(modelLandmarks.LM,modelLandmarks.MM);
modelLandmarks.IC = midPoint3d(modelLandmarks.LC,modelLandmarks.MC);

%Create variable for model landmark names
modelLandmarkNames = fieldnames(modelLandmarks);

%Create new landmarks on tibial anterior and fibula lateral border
[modelLandmarks, modelLandmarkNames, modelLandmarkInds] = createBoneBorderLandmarks_shapeModel(...
    tibiaFibulaShapeModel, modelLandmarks, modelLandmarkNames, palpableLandmarks, ...
    optionStruct_tib.nb_pts, optionStruct_fib.nb_pts);

% % % %Visualise mean shape model and landmarks
% % % cFigure; hold on
% % % gpatch(tibiaFibulaShapeModel.F,tibiaFibulaShapeModel.meanPoints,'gw','k', 0.5)
% % % for landmarkNo = 1:length(modelLandmarkNames)
% % %     plotV(modelLandmarks.(modelLandmarkNames{landmarkNo}), 'r.','MarkerSize', 25);
% % % end
% % % axisGeom; camlight headlight

%% Load and reconstruct sample data

%Get sample ID's from folder
sampleDir = dir('samples\*.stl');
for sampleId = 1:length(sampleDir)
    sampleNames{sampleId} = sampleDir(sampleId).name;
end

%Check whether to run reconstructions
if runReconstructions

    %Set-up waitbar for processing surfaces
    wbar_surfaces = waitbar(0, 'Reconstructing surfaces from landmarks...');
    
    %Loop through samples
    for sampleNo = 1:length(sampleNames)

        %Set participant Id
        pID = strsplit(sampleNames{sampleNo},'.');
        pID = pID{1};

        %Tibia-fibula
        [tibFibSTLstruct] = import_STL(['samples\',pID,'.stl']);
        tibFibF = tibFibSTLstruct.solidFaces{1}; %Faces
        tibFibV = tibFibSTLstruct.solidVertices{1}; %Vertices
        [tibFibF,tibFibV] = mergeVertices(tibFibF,tibFibV);

        %Split the surfaces to grab the tibia only

        %Group the vertices and faces from the combined shape model
        [groupIndexVertices,groupIndexFaces] = groupVertices(tibFibF, tibFibV,0);

        %Identify which grouped section contains a higher volume
        %This will be indicative of the tibia
        %%%%% NOTE: this works for every sample except C32RTF for seemingly
        %%%%% the fact that the tibia is so thin...
        if tetVolMeanEst(tibFibF(groupIndexFaces == 1,:),tibFibV) > ...
                tetVolMeanEst(tibFibF(groupIndexFaces == 2,:),tibFibV)
            %First index is tibia
            logicKeep = groupIndexFaces == 1;
        else
            %Second index is tibia
            logicKeep = groupIndexFaces == 2;
        end
        %Separate the surfaces
        if sampleNo ~= 32
            [tibiaF, tibiaV] = patchCleanUnused(tibFibF(logicKeep,:), tibFibV);
            [fibulaF, fibulaV] = patchCleanUnused(tibFibF(~logicKeep,:), tibFibV);
        else
            [tibiaF, tibiaV] = patchCleanUnused(tibFibF(~logicKeep,:), tibFibV);
            [fibulaF, fibulaV] = patchCleanUnused(tibFibF(logicKeep,:), tibFibV);
        end

% % %         %Visualise
% % %         cFigure; hold on
% % %         gpatch(tibiaF,tibiaV,'gw', 'none')
% % %         gpatch(fibulaF,fibulaV,'bw', 'none')
% % %         axisGeom; camlight headlight

        %Remesh surfaces

        %Run initial smoothing on surfaces
        smoothPar.n = 2; %2 iterations
        smoothPar.Method = 'LAP'; %Laplacian smoothing method        
        tibiaV = patchSmooth(tibiaF,tibiaV,[],smoothPar);
        fibulaV = patchSmooth(fibulaF,fibulaV,[],smoothPar);

        %Remesh using ggremesh to reduce complexity
        [tibiaF,tibiaV] = ggremesh(tibiaF, tibiaV, optionStruct_tib);

        %Check for matching point number and display if not
        if length(tibiaV) ~= optionStruct_tib.nb_pts
            error('Number of tibia vertices does not match requested for case %s\n',caseNo);
        end

        %Remesh using ggremesh to reduce complexity
        [fibulaF,fibulaV] = ggremesh(fibulaF,fibulaV, optionStruct_fib);

        %Check for matching point number and display if not
        if length(fibulaV) ~= optionStruct_fib.nb_pts
            %Print out summary this time
            error('Number of fibula vertices does not match requested for case %s\n',caseNo);
        end

        %Check for holes in surfaces
        %This shouldn't occur with the remeshing, but check for safety...
        if ~isempty(patchBoundary(tibiaF, tibiaV))
            error('Holes in tibia mesh')
        end
        if ~isempty(patchBoundary(fibulaF, fibulaV))
            error('Holes in fibula mesh')
        end

        %Load sample landmarks
        sampleLandmarkData = readtable(['samples\',pID,'.csv']);
        sampleLandmarkData.Properties.VariableNames = {'landmark' 'X' 'Y' 'Z'};

        %Extract landmarks into structure as variables
        for landmarkNo = 1:length(sampleLandmarkData.landmark)
            sampleLandmarks.(sampleLandmarkData.landmark{landmarkNo}) = ...
                [sampleLandmarkData.X(landmarkNo), ...
                sampleLandmarkData.Y(landmarkNo), ...
                sampleLandmarkData.Z(landmarkNo)];    
        end

        %Create new landmarks for mean model
        sampleLandmarks.IM = midPoint3d(sampleLandmarks.LM,sampleLandmarks.MM);
        sampleLandmarks.IC = midPoint3d(sampleLandmarks.LC,sampleLandmarks.MC);

        %Create variable for model landmark names
        sampleLandmarkNames = fieldnames(sampleLandmarks);

        %Rotate the landmarks to be in the tibial coordinate system like the shape model
        [sampleLandmarks, tibiaV, fibulaV] = ...
            alignLandmarksSurfaces(sampleLandmarks, tibiaV, fibulaV);

% % %         %Visualise sample model and landmarks
% % %         cFigure; hold on
% % %         gpatch(tibiaF,tibiaV,'gw','k', 0.5)
% % %         gpatch(fibulaF,fibulaV,'bw','k', 0.5)
% % %         for landmarkNo = 1:length(sampleLandmarkNames)
% % %             plotV(sampleLandmarks.(sampleLandmarkNames{landmarkNo}), 'r.','MarkerSize', 25);
% % %         end
% % %         axisGeom; camlight headlight

        %Create new landmarks on tibial anterior and fibula lateral border for sample
        [sampleLandmarks, sampleLandmarkNames] = createBoneBorderLandmarks_sample(...
            sampleLandmarks, sampleLandmarkNames, palpableLandmarks, ...
            tibiaV, fibulaV);

% % %         %Visualise rotated surfaces and landmarks
% % %         cFigure; hold on;
% % %         %Tibia
% % %         gpatch(tibiaF,tibiaV,'gw', 'k')
% % %         %Fibula
% % %         gpatch(fibulaF,fibulaV,'bw', 'k')
% % %         %Palpable landmarks
% % %         for landmarkNo = 1:length(palpableLandmarks)
% % %             plotV(sampleLandmarks.(palpableLandmarks{landmarkNo}), 'r.', 'MarkerSize', 25)
% % %         end
% % %         %Axis parameters
% % %         axisGeom; camlight headlight

        %% Optimise the shape model parameters to fit the dataset

        %Create the function handle for the optimisation
        optFunc = @(pcScores)calcLandmarkError(pcScores, tibiaFibulaShapeModel, sampleLandmarks, modelLandmarkInds);

% % %         %Register the new surface against the shape model mean to create an
% % %         %appropriate initial guess
% % % 
% % %         %Rigidly register the tibia and fibula to the mean of the shape model
% % %         %to perform rough initial alignment
% % %         [cpdTformRig] = cpd_register(tibiaFibulaShapeModel.meanPoints, [tibiaV; fibulaV], optCPD_rig);
% % % 
% % %         %Replace the tibia and fibula points with the rigidly registered
% % %         tibiaV = cpdTformRig.Y(1:optionStruct_tib.nb_pts,:);
% % %         fibulaV = cpdTformRig.Y(optionStruct_tib.nb_pts+1:end,:);
% % % 
% % %         %Perform non-rigid registration using CPD to map new surface to shape model
% % %         [cpdTformIG_tib] = cpd_register(tibiaFibulaShapeModel.meanPoints(1:optionStruct_tib.nb_pts,:), tibiaV, optCPD);
% % %         [cpdTformIG_fib] = cpd_register(tibiaFibulaShapeModel.meanPoints(optionStruct_tib.nb_pts+1:end,:), fibulaV, optCPD);
% % % 
% % %         %Identify the matching points in the target mesh against the
% % %         %registration to identify the corresponding point indices
% % %         regSortIdxIG_tib = knnsearch(cpdTformIG_tib.Y, tibiaFibulaShapeModel.meanPoints(1:optionStruct_tib.nb_pts,:));
% % %         regSortIdxIG_fib = knnsearch(cpdTformIG_fib.Y, tibiaFibulaShapeModel.meanPoints(optionStruct_tib.nb_pts+1:end,:));
% % % 
% % %         %Sort the registered points so they align with the reference mesh
% % %         igDataReg_tib = tibiaV(regSortIdxIG_tib,:);
% % %         igDataReg_fib = fibulaV(regSortIdxIG_fib,:);
% % % 
% % %         %Reshape the registered data points and remove the mean
% % %         igData = reshape([igDataReg_tib;igDataReg_fib]', ...
% % %             [1, length([igDataReg_tib;igDataReg_fib])*3]) - tibiaFibulaShapeModel.mean;
% % %         
% % %         %Project the new data against the loadings to find the estimated PC
% % %         %scores for each component
% % %         for nPC = 1:tibiaFibulaShapeModel.retainPCs
% % %             %Calculate scores and set to an initial guess variable
% % %             x0(nPC,1) = dot(tibiaFibulaShapeModel.loadings(:,nPC), igData);
% % %         end

% % %         %Project the new data against the loadings to find the estimated scores
% % %         %Set the upper and lower bounds on the PC score variables within the optimisation
% % %         %Here we set them as -5 / +5 standard deviations about the PC score mean (i.e. zero)
% % %         %An initial guess of the mean (i.e. zero) is also set here
% % %         sdRange = 5;
% % %         %Loop through retained PCs
% % %         for nPC = 1:tibiaFibulaShapeModel.retainPCs
% % %             %Calculate scores and set to an initial guess variable
% % %             x0(nPC,1) = dot(tibiaFibulaShapeModel.loadings(:,nPC), igData);
% % %             %Set the lower bound
% % %             lb(nPC,1) = std(tibiaFibulaShapeModel.score(:,nPC)) * -sdRange;
% % %             %Set the upper bound
% % %             ub(nPC,1) = std(tibiaFibulaShapeModel.score(:,nPC)) * sdRange;
% % %         end

        %Set the upper and lower bounds on the PC score variables within the optimisation
        %Here we set them as -5 / +5 standard deviations about the PC score mean (i.e. zero)
        %An initial guess of the mean (i.e. zero) is also set here
        sdRange = 5;
        %Loop through retained PCs
        for nPC = 1:tibiaFibulaShapeModel.retainPCs
            %Calculate scores and set to an initial guess variable
            x0(nPC,1) = 0;
            %Set the lower bound
            lb(nPC,1) = std(tibiaFibulaShapeModel.score(:,nPC)) * -sdRange;
            %Set the upper bound
            ub(nPC,1) = std(tibiaFibulaShapeModel.score(:,nPC)) * sdRange;
        end

        %Set fmincon options
        fminconOpts = optimoptions('fmincon');
        fminconOpts.MaxIterations = 5000; %up potential iterations
        fminconOpts.MaxFunctionEvaluations = 20000; %up potential function evals
        fminconOpts.StepTolerance = 1e-4; %scale optimisation tolerance to problem
        fminconOpts.Display = 'iter'; %set to display function iterations

        %Run the optimisation to get the PC scores that miniises the difference
        %between the sample and shape model landmarks
        [optPcScores, landmarkError] = fmincon(optFunc, x0, [], [], [], [], lb, ub, [], fminconOpts);

% % %         %Reconstruct the surfaces using the projected PC scores
% % %         optReconstructed = x0(1:tibiaFibulaShapeModel.retainPCs)' * ...
% % %             tibiaFibulaShapeModel.loadings(:,1:tibiaFibulaShapeModel.retainPCs)' + ...
% % %             tibiaFibulaShapeModel.mean;

        %Reconstruct using the optimised PC scores
        optReconstructed = optPcScores(1:tibiaFibulaShapeModel.retainPCs)' * ...
            tibiaFibulaShapeModel.loadings(:,1:tibiaFibulaShapeModel.retainPCs)' + ...
            tibiaFibulaShapeModel.mean;

        %Reshape reconstructed points
        optReconstructedV = reshape(optReconstructed', [3, length(optReconstructed)/3])';
        
        %Extract the transformation used to align the landmarks
        %This can then be applied to the surface to calculate error
        %Extract shape model and reconstruct landmark locations used in the reconstruction
        reconLandmarkNames = fieldnames(modelLandmarkInds);
        for landmarkNo = 1:length(reconLandmarkNames)
            reconLandmarks(landmarkNo,:) = optReconstructedV(modelLandmarkInds.(reconLandmarkNames{landmarkNo}),:);
            origLandmarks(landmarkNo,:) = sampleLandmarks.(reconLandmarkNames{landmarkNo});
        end
        %Get the transformation
        [~,reconAlignedLandmarks,tform] = procrustes(origLandmarks, reconLandmarks,...
            'Scaling', false, 'Reflection', false);
        
        %Align the reconstructed points using the transformation
        optReconstructedV = tform.b * optReconstructedV * tform.T + mean(tform.c);

        %Unmerge the tibia and fibula in the reconstruction
        tibiaReconstructedV = optReconstructedV(1:optionStruct_tib.nb_pts,:);
        fibulaReconstructedV = optReconstructedV(optionStruct_tib.nb_pts+1:end,:);

        %Register the reconstructions to the original surface
        [cpdTformTib,regSortIdxTib] = cpd_register(tibiaReconstructedV, tibiaV, optCPD);
        [cpdTformFib,regSortIdxFib] = cpd_register(fibulaReconstructedV, fibulaV, optCPD);

% % %         %Identify the matching points in the target mesh against the
% % %         %registration to identify the corresponding point indices
% % %         regSortIdxTib = knnsearch(cpdTformTib.Y, tibiaV);
% % %         regSortIdxFib = knnsearch(cpdTformFib.Y, fibulaV);

        %Sort the registered points so they align with the reference mesh
        tibiaReconstructedV = tibiaReconstructedV(regSortIdxTib,:);
        fibulaReconstructedV = fibulaReconstructedV(regSortIdxFib,:);
        
% % %         %Rigidly register the reconstructed tibia and fibula to the
% % %         %original surface to quantify error
% % %         [cpdTformRig] = cpd_register([tibiaV; fibulaV], ...
% % %             [tibiaReconstructedV; fibulaReconstructedV], optCPD_rig);
% % %         
% % %         %Replace the tibia and fibula points with the rigidly registered
% % %         tibiaReconstructedV = cpdTformRig.Y(1:optionStruct_tib.nb_pts,:);
% % %         fibulaReconstructedV = cpdTformRig.Y(optionStruct_tib.nb_pts+1:end,:);

% % %         %Visualise original vs. reconstructed
% % %         cFigure; hold on;
% % %         %Originals
% % %         gpatch(tibiaF, tibiaV, 'gw', 'none', 0.3);
% % %         gpatch(fibulaF, fibulaV, 'gw', 'none', 0.3);
% % %         %Reconstructions
% % %         gpatch(tibiaF, tibiaReconstructedV, 'rw', 'k', 1);
% % %         gpatch(fibulaF, fibulaReconstructedV, 'rw', 'k', 1);
% % %         axisGeom; camlight headlight

        %% Quantify and visualise the error of the reconstruction
        
        %Add sum and mean of landmark error
        landmarkErrorDistSum(sampleNo,1) = sum(distancePoints3d(reconAlignedLandmarks, origLandmarks));
        landmarkErrorDistMean(sampleNo,1) = mean(distancePoints3d(reconAlignedLandmarks, origLandmarks));

        %Calculate point error distance
        pointErrorDistTib(sampleNo,:) = distancePoints3d(tibiaV, tibiaReconstructedV);
        pointErrorDistFib(sampleNo,:) = distancePoints3d(fibulaV, fibulaReconstructedV);

        %Calculate the mean error
        pointErrorDistTibMean(sampleNo,1) = mean(pointErrorDistTib(sampleNo,:));
        pointErrorDistFibMean(sampleNo,1) = mean(pointErrorDistFib(sampleNo,:));

        %Calculate the peak error
        pointErrorDistTibMax(sampleNo,1) = max(pointErrorDistTib(sampleNo,:));
        pointErrorDistFibMax(sampleNo,1) = max(pointErrorDistFib(sampleNo,:));

        %Convert distance error to colour scales for visualisation
        pointErrorDistTibColF(sampleNo,:) = vertexToFaceMeasure(tibiaF, pointErrorDistTib(sampleNo,:)');
        pointErrorDistTibColV(sampleNo,:) = faceToVertexMeasure(tibiaF, tibiaV, pointErrorDistTibColF(sampleNo,:)');
        pointErrorDistFibColF(sampleNo,:) = vertexToFaceMeasure(fibulaF, pointErrorDistFib(sampleNo,:)');
        pointErrorDistFibColV(sampleNo,:) = faceToVertexMeasure(fibulaF, fibulaV, pointErrorDistFibColF(sampleNo,:)');

        %Calculate Jaccard similarity across both parts and overall
        [jaccardSimilarityTibia(sampleNo,1), jaccardSimilarityFibula(sampleNo,1), jaccardSimilarityAll(sampleNo,1)] = ...
            calcJaccardLandmarkReconstructed(tibiaF, fibulaF, tibiaV, fibulaV, tibiaReconstructedV, fibulaReconstructedV);

        %Create visualisation of error
        %Use subplots to create different perspectives
        cFigure; hold on;
        %Loop through four views to create subplot
        for viewNo = 1:4    
            %Create subplot for current view
            subplot(1,4,viewNo);
            %Add greyed out original surface
            hpOrigTib = gpatch(tibiaF, tibiaV, [200/255 200/255 200/255], 'none', 0.5);
            hpOrigFib = gpatch(fibulaF, fibulaV, [200/255 200/255 200/255], 'none', 0.5);
            %Add colormapped reconstruction
            hpPredTib = gpatch(tibiaF, tibiaReconstructedV, pointErrorDistTibColV(sampleNo,:)', 'none', 1);
            hpPredFib = gpatch(fibulaF, fibulaReconstructedV, pointErrorDistFibColV(sampleNo,:)', 'none', 1);
            %Interpolate colouring for smoothness
            hpPredTib.FaceColor = 'Interp';
            hpPredFib.FaceColor = 'Interp';
            colormap viridis
            %Set axis view
            axis equal; axis tight; view(0,90);
            rotate(hpOrigTib,[0 1 0], surfaceRot(viewNo));
            rotate(hpOrigFib,[0 1 0], surfaceRot(viewNo));
            rotate(hpPredTib,[0 1 0], surfaceRot(viewNo));
            rotate(hpPredFib,[0 1 0], surfaceRot(viewNo));
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
        export_fig(['figures\reconstructionError\',pID,'_reconstructionErrorMap.png'],'-m1');
        close

        %Save predicted surface as an STL
        %Tibia
        stlStructTib.solidNames = {[pID,'-predicted-tibia']}; %names of parts
        stlStructTib.solidVertices = {tibiaReconstructedV}; %Vertices
        stlStructTib.solidFaces = {tibiaF}; %Faces
        stlStructTib.solidNormals={[]};
        %Fibula
        stlStructFib.solidNames = {[pID,'-predicted-fibula']}; %names of parts
        stlStructFib.solidVertices = {fibulaReconstructedV}; %Vertices
        stlStructFib.solidFaces = {fibulaF}; %Faces
        stlStructFib.solidNormals={[]};

        %Export STLs
        export_STL_txt(['predictedSurfaces\',pID,'-predicted-tibia.stl'], stlStructTib);
        export_STL_txt(['predictedSurfaces\',pID,'-predicted-fibula.stl'], stlStructFib);

        %Store data in array
        predictedSurfaces.tibiaV(:,:,sampleNo) = tibiaReconstructedV;
        predictedSurfaces.tibiaF(:,:,sampleNo) = tibiaF;
        predictedSurfaces.fibulaV(:,:,sampleNo) = fibulaReconstructedV;
        predictedSurfaces.fibulaF(:,:,sampleNo) = fibulaF;

        %Clear out the landmarks as these don't get reallocated
        clear sampleLandmarks

        %Update waitbar
        waitbar(sampleNo/length(sampleNames), wbar_surfaces);

    %% End of sample loop
    end

    %Close waitbar
    close(wbar_surfaces);

    %Place calculations into structure and tables
    reconstructionErrorSummary.errorSummaryTable = table(landmarkErrorDistSum, ...
        landmarkErrorDistMean, ...
        pointErrorDistTibMean, ...
        pointErrorDistFibMean, ...
        pointErrorDistTibMax, ...
        pointErrorDistFibMax, ...
        jaccardSimilarityTibia, ...
        jaccardSimilarityFibula, ...
        jaccardSimilarityAll, ...
        'VariableNames', {'landmarkErrorDistSum', 'landmarkErrorDistMean', 'pointErrorDistTibMean', 'pointErrorDistFibMean', 'pointErrorDistTibMax', ...
        'pointErrorDistFibMax', 'jaccardSimilarityTibia', 'jaccardSimilarityFibula', 'jaccardSimilarityAll'}, ...
        'RowNames', sampleNames);
    
    %%%%% TODO: review the second run through of results...working better
    %%%%% without the optimisation, but still some issues mainly with
    %%%%% fibula volumetric meshing to calculate Jaccard value

    %Export table to file
    writetable(reconstructionErrorSummary.errorSummaryTable, 'results\errorSummaryTable.csv', ...
        'WriteRowNames', true);

    %Place results in structure
    reconstructionErrorSummary.landmarkErrorDistSum = landmarkErrorDistSum;
    reconstructionErrorSummary.landmarkErrorDistMean = landmarkErrorDistMean;
    reconstructionErrorSummary.jaccardSimilarityTibia = jaccardSimilarityTibia;
    reconstructionErrorSummary.jaccardSimilarityFibula = jaccardSimilarityFibula;
    reconstructionErrorSummary.jaccardSimilarityAll = jaccardSimilarityAll;
    reconstructionErrorSummary.pointErrorDistTib = pointErrorDistTib;
    reconstructionErrorSummary.pointErrorDistFib = pointErrorDistFib;
    reconstructionErrorSummary.pointErrorDistTibMean = pointErrorDistTibMean;
    reconstructionErrorSummary.pointErrorDistFibMean = pointErrorDistFibMean;
    reconstructionErrorSummary.pointErrorDistTibMax = pointErrorDistTibMax;
    reconstructionErrorSummary.pointErrorDistFibMax = pointErrorDistFibMax;
    reconstructionErrorSummary.pointErrorDistTibColF = pointErrorDistTibColF;
    reconstructionErrorSummary.pointErrorDistFibColF = pointErrorDistFibColF;
    reconstructionErrorSummary.pointErrorDistTibColV = pointErrorDistTibColV;
    reconstructionErrorSummary.pointErrorDistFibColV = pointErrorDistFibColV;
    reconstructionErrorSummary.predictedSurfaces.tibiaV = predictedSurfaces.tibiaV;
    reconstructionErrorSummary.predictedSurfaces.fibulaV = predictedSurfaces.fibulaV;
    reconstructionErrorSummary.predictedSurfaces.tibiaF = predictedSurfaces.tibiaF;
    reconstructionErrorSummary.predictedSurfaces.fibulaF = predictedSurfaces.fibulaF;

    %Save error summary as mat
    save('results\reconstructionErrorSummary.mat', 'reconstructionErrorSummary');

else
    
    %Load pre-processed results
    load('results\reconstructionErrorSummary.mat');
    
end

%Create summary figure of error data alongside calculations
%Create figure
hfErr = cFigure; hold on
hfErr.Units = 'centimeters';
hfErr.Position = [5, 5, 24, 18];

%Calculate and display mean and 95% CI's for the mean and peak error

%All
%Mean & SD for mean landmark error
meanError_m = mean(reconstructionErrorSummary.landmarkErrorDistMean);
meanError_sd = std(reconstructionErrorSummary.landmarkErrorDistMean);
meanError_lower95 = meanError_m - (1.96 * (meanError_sd / sqrt(length(sampleNames))));
meanError_upper95 = meanError_m + (1.96 * (meanError_sd / sqrt(length(sampleNames))));
%Mean & SD for landmark error sum
sumError_m = mean(reconstructionErrorSummary.landmarkErrorDistSum);
sumError_sd = std(reconstructionErrorSummary.landmarkErrorDistSum);
sumError_lower95 = sumError_m - (1.96 * (sumError_sd / sqrt(length(sampleNames))));
sumError_upper95 = sumError_m + (1.96 * (sumError_sd / sqrt(length(sampleNames))));
%Jaccard similarity
jaccard_m = mean(reconstructionErrorSummary.jaccardSimilarityAll);
jaccard_sd = std(reconstructionErrorSummary.jaccardSimilarityAll);
jaccard_lower95 = jaccard_m - (1.96 * (jaccard_sd / sqrt(length(sampleNames))));
jaccard_upper95 = jaccard_m + (1.96 * (jaccard_sd / sqrt(length(sampleNames))));
%Display
disp(['Mean error in landmark point distance (mean +/- 95% CIs) = ',num2str(round(meanError_m,2)), ...
    '[',num2str(round(meanError_lower95,2)),',',num2str(round(meanError_upper95,2)),']']);
disp(['Sum error in landmark point distance (mean +/- 95% CIs) = ',num2str(round(sumError_m,2)), ...
    '[',num2str(round(sumError_lower95,2)),',',num2str(round(sumError_upper95,2)),']']);
disp(['Overall reconstruction Jaccard Index (mean +/- 95% CIs) = ',num2str(round(jaccard_m,3)), ...
    '[',num2str(round(jaccard_lower95,3)),',',num2str(round(jaccard_upper95,3)),']']);

%Plot mean landmark error data
%Create subplot
subplot(3,3,1); hold on
%Plot mean as dashed line & CIs as dotted lines
plot([1,length(sampleNames)], [meanError_m, meanError_m], ...
    'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
plot([1,length(sampleNames)], [meanError_lower95, meanError_lower95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
plot([1,length(sampleNames)], [meanError_upper95, meanError_upper95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
%Plot participant data
scatter(linspace(1,length(sampleNames), length(sampleNames)), reconstructionErrorSummary.landmarkErrorDistMean, ...
    15, 'k', 'filled');
%Set labels
xlabel('Samples'); ylabel('Error (mm)'); title('Mean Landmark Error');
%Remove ticks
ax = gca(); ax.XTick = []; ax.XLim = [0,37];
%Plot sum landmark error data
%Create subplot
subplot(3,3,2); hold on
%Plot mean as dashed line & CIs as dotted lines
plot([1,length(sampleNames)], [sumError_m, sumError_m], ...
    'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
plot([1,length(sampleNames)], [sumError_lower95, sumError_lower95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
plot([1,length(sampleNames)], [sumError_upper95, sumError_upper95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
%Plot participant data
scatter(linspace(1,length(sampleNames), length(sampleNames)), reconstructionErrorSummary.landmarkErrorDistSum, ...
    15, 'k', 'filled');
%Set labels
xlabel('Samples'); ylabel('Error (mm)'); title('Sum Landmark Error');
%Remove ticks
ax = gca(); ax.XTick = []; ax.XLim = [0,37];
%Plot Jaccard similarity data
%Create subplot
subplot(3,3,3); hold on
%Plot mean as dashed line & CIs as dotted lines
plot([1,length(sampleNames)], [jaccard_m, jaccard_m], ...
    'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
plot([1,length(sampleNames)], [jaccard_lower95, jaccard_lower95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
plot([1,length(sampleNames)], [jaccard_upper95, jaccard_upper95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
%Plot participant data
scatter(linspace(1,length(sampleNames), length(sampleNames)), reconstructionErrorSummary.jaccardSimilarityAll, ...
    15, 'k', 'filled');
%Set labels
xlabel('Samples'); ylabel('Jaccard Index (0-1)'); title('Jaccard Similarity (All)');
%Remove ticks
ax = gca(); ax.XTick = []; ax.XLim = [0,37];

%Tibia
%Mean & SD for mean error
meanError_m = mean(reconstructionErrorSummary.pointErrorDistTibMean);
meanError_sd = std(reconstructionErrorSummary.pointErrorDistTibMean);
meanError_lower95 = meanError_m - (1.96 * (meanError_sd / sqrt(length(sampleNames))));
meanError_upper95 = meanError_m + (1.96 * (meanError_sd / sqrt(length(sampleNames))));
%Mean & SD for maxerror
maxError_m = mean(reconstructionErrorSummary.pointErrorDistTibMax);
maxError_sd = std(reconstructionErrorSummary.pointErrorDistTibMax);
maxError_lower95 = maxError_m - (1.96 * (maxError_sd / sqrt(length(sampleNames))));
maxError_upper95 = maxError_m + (1.96 * (maxError_sd / sqrt(length(sampleNames))));
%Jaccard similarity
jaccard_m = mean(reconstructionErrorSummary.jaccardSimilarityTibia);
jaccard_sd = std(reconstructionErrorSummary.jaccardSimilarityTibia);
jaccard_lower95 = jaccard_m - (1.96 * (jaccard_sd / sqrt(length(sampleNames))));
jaccard_upper95 = jaccard_m + (1.96 * (jaccard_sd / sqrt(length(sampleNames))));
%Display
disp(['Tibia reconstruction mean point error (mean +/- 95% CIs) = ',num2str(round(meanError_m,2)), ...
    '[',num2str(round(meanError_lower95,2)),',',num2str(round(meanError_upper95,2)),']']);
disp(['Tibia reconstruction max point error (mean +/- 95% CIs) = ',num2str(round(maxError_m,2)), ...
    '[',num2str(round(maxError_lower95,2)),',',num2str(round(maxError_upper95,2)),']']);
disp(['Tibia reconstruction Jaccard Index (mean +/- 95% CIs) = ',num2str(round(jaccard_m,3)), ...
    '[',num2str(round(jaccard_lower95,3)),',',num2str(round(jaccard_upper95,3)),']']);

%Plot mean point error
%Create subplot
subplot(3,3,4); hold on
%Plot mean as dashed line & CIs as dotted lines
plot([1,length(sampleNames)], [meanError_m, meanError_m], ...
    'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
plot([1,length(sampleNames)], [meanError_lower95, meanError_lower95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
plot([1,length(sampleNames)], [meanError_upper95, meanError_upper95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
%Plot participant data
scatter(linspace(1,length(sampleNames), length(sampleNames)), reconstructionErrorSummary.pointErrorDistTibMean, ...
    15, 'k', 'filled');
%Set labels
xlabel('Samples'); ylabel('Error (mm)'); title('Mean Position Error (Tibia)');
%Remove ticks
ax = gca(); ax.XTick = []; ax.XLim = [0,37];
%Plot max point error
%Create subplot
subplot(3,3,5); hold on
%Plot mean as dashed line & CIs as dotted lines
plot([1,length(sampleNames)], [maxError_m, maxError_m], ...
    'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
plot([1,length(sampleNames)], [maxError_lower95, maxError_lower95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
plot([1,length(sampleNames)], [maxError_upper95, maxError_upper95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
%Plot participant data
scatter(linspace(1,length(sampleNames), length(sampleNames)), reconstructionErrorSummary.pointErrorDistTibMax, ...
    15, 'k', 'filled');
%Set labels
xlabel('Samples'); ylabel('Error (mm)'); title('Peak Position Error (Tibia)');
%Remove ticks
ax = gca(); ax.XTick = []; ax.XLim = [0,37];
%Plot Jaccard similarity data
%Create subplot
subplot(3,3,6); hold on
%Plot mean as dashed line & CIs as dotted lines
plot([1,length(sampleNames)], [jaccard_m, jaccard_m], ...
    'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
plot([1,length(sampleNames)], [jaccard_lower95, jaccard_lower95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
plot([1,length(sampleNames)], [jaccard_upper95, jaccard_upper95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
%Plot participant data
scatter(linspace(1,length(sampleNames), length(sampleNames)), reconstructionErrorSummary.jaccardSimilarityTibia, ...
    15, 'k', 'filled');
%Set labels
xlabel('Samples'); ylabel('Jaccard Index (0-1)'); title('Jaccard Similarity (Tibia)');
%Remove ticks
ax = gca(); ax.XTick = []; ax.XLim = [0,37];

%Fibula
%Mean & SD for mean error
meanError_m = mean(reconstructionErrorSummary.pointErrorDistFibMean);
meanError_sd = std(reconstructionErrorSummary.pointErrorDistFibMean);
meanError_lower95 = meanError_m - (1.96 * (meanError_sd / sqrt(length(sampleNames))));
meanError_upper95 = meanError_m + (1.96 * (meanError_sd / sqrt(length(sampleNames))));
%Mean & SD for maxerror
maxError_m = mean(reconstructionErrorSummary.pointErrorDistFibMax);
maxError_sd = std(reconstructionErrorSummary.pointErrorDistFibMax);
maxError_lower95 = maxError_m - (1.96 * (maxError_sd / sqrt(length(sampleNames))));
maxError_upper95 = maxError_m + (1.96 * (maxError_sd / sqrt(length(sampleNames))));
%Jaccard similarity
jaccard_m = mean(reconstructionErrorSummary.jaccardSimilarityFibula);
jaccard_sd = std(reconstructionErrorSummary.jaccardSimilarityFibula);
jaccard_lower95 = jaccard_m - (1.96 * (jaccard_sd / sqrt(length(sampleNames))));
jaccard_upper95 = jaccard_m + (1.96 * (jaccard_sd / sqrt(length(sampleNames))));
%Display
disp(['Fibula reconstruction mean point error (mean +/- 95% CIs) = ',num2str(round(meanError_m,2)), ...
    '[',num2str(round(meanError_lower95,2)),',',num2str(round(meanError_upper95,2)),']']);
disp(['Fibula reconstruction max point error (mean +/- 95% CIs) = ',num2str(round(maxError_m,2)), ...
    '[',num2str(round(maxError_lower95,2)),',',num2str(round(maxError_upper95,2)),']']);
disp(['Fibula reconstruction Jaccard Index (mean +/- 95% CIs) = ',num2str(round(jaccard_m,3)), ...
    '[',num2str(round(jaccard_lower95,3)),',',num2str(round(jaccard_upper95,3)),']']);

%Plot mean point error
%Create subplot
subplot(3,3,7); hold on
%Plot mean as dashed line & CIs as dotted lines
plot([1,length(sampleNames)], [meanError_m, meanError_m], ...
    'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
plot([1,length(sampleNames)], [meanError_lower95, meanError_lower95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
plot([1,length(sampleNames)], [meanError_upper95, meanError_upper95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
%Plot participant data
scatter(linspace(1,length(sampleNames), length(sampleNames)), reconstructionErrorSummary.pointErrorDistFibMean, ...
    15, 'k', 'filled');
%Set labels
xlabel('Samples'); ylabel('Error (mm)'); title('Mean Position Error (Fibula)');
%Remove ticks
ax = gca(); ax.XTick = []; ax.XLim = [0,37];
%Plot max point error
%Create subplot
subplot(3,3,8); hold on
%Plot mean as dashed line & CIs as dotted lines
plot([1,length(sampleNames)], [maxError_m, maxError_m], ...
    'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
plot([1,length(sampleNames)], [maxError_lower95, maxError_lower95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
plot([1,length(sampleNames)], [maxError_upper95, maxError_upper95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
%Plot participant data
scatter(linspace(1,length(sampleNames), length(sampleNames)), reconstructionErrorSummary.pointErrorDistFibMax, ...
    15, 'k', 'filled');
%Set labels
xlabel('Samples'); ylabel('Error (mm)'); title('Peak Position Error (Fibula)');
%Remove ticks
ax = gca(); ax.XTick = []; ax.XLim = [0,37];
%Plot Jaccard similarity data
%Create subplot
subplot(3,3,9); hold on
%Plot mean as dashed line & CIs as dotted lines
plot([1,length(sampleNames)], [jaccard_m, jaccard_m], ...
    'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
plot([1,length(sampleNames)], [jaccard_lower95, jaccard_lower95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
plot([1,length(sampleNames)], [jaccard_upper95, jaccard_upper95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
%Plot participant data
scatter(linspace(1,length(sampleNames), length(sampleNames)), reconstructionErrorSummary.jaccardSimilarityFibula, ...
    15, 'k', 'filled');
%Set labels
xlabel('Samples'); ylabel('Jaccard Index (0-1)'); title('Jaccard Similarity (Fibula)');
%Remove ticks
ax = gca(); ax.XTick = []; ax.XLim = [0,37];

%Export summary figure
export_fig('figures\errorSummary\errorSummaryFigure.png','-m2');
close(hfErr);

%% ----- end of generateSurfaceFromLandmarks.m ----- %%