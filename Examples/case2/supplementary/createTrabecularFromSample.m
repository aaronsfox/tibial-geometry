function [newPredictedTrabV, newTibiaF, newTibiaV] = createTrabecularFromSample(tibiaF, tibiaV, tibTrabShapeModel, nPtsTibia)

    %% This script creates a trabecular surface from one of the sample
    %  tibia-fibula surfaces from the Nolte et al. dataset
    %
    %  Inputs:
    %       tibiaF - faces for the sample tibia surface
    %       tibiaV - vertices for the sample tibia surface
    %       tibTrabShapeModel - the structure with shape model data for the tibia + trabecular
    %       nPtsTibia - number of vertices on tibial surface in shape model
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

    %% Set-up

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
    
    %% Generate new trabecular surface

    %It helps the registration algorithm to get the surfaces roughly aligned, 
    %so perform some initial transformations to get these objects in a similar
    %orientation and position. The rotations and transformations applied here 
    %were manually identified through some eye-balling.

    %Rotate -90 degrees about the x-axis
    tibiaV = transformPoint3d(tibiaV, createRotationOx(deg2rad(-90)));  

    %Rotate 90 degrees about the y-axis
    tibiaV = transformPoint3d(tibiaV, createRotationOy(deg2rad(90)));

    %Identify the magnitude required for a rough translation by comparing the
    %means of the points
    %Create translation matrix
    transMat = createTranslation3d(mean(tibTrabShapeModel.meanPoints(1:nPtsTibia,:)) - mean(tibiaV));
    %Apply translation to points
    tibiaV = transformPoint3d(tibiaV,transMat);

    %Apply a rigid transformation to match the new tibia points with the shape model
    %Don't apply the scaling from this rigid transform
    [cpdTformTibRig] = cpd_register(tibTrabShapeModel.meanPoints(1:nPtsTibia,:), tibiaV, optCPD_rig);
    tibiaV = tibiaV * cpdTformTibRig.R' + repmat(cpdTformTibRig.t,1,length(tibiaV))';

    %To map the new surface to the shape model it's important to identify the
    %corresponding points on the surface to that of the shape model.

    %Perform non-rigid registration using CPD to align surfaces
    [cpdTformTib] = cpd_register(tibTrabShapeModel.meanPoints(1:nPtsTibia,:), tibiaV, optCPD);
    
    %Identify the matching points in the target mesh against the
    %registration to identify the corresponding point indices
    regSortIdxTib = knnsearch(cpdTformTib.Y, tibTrabShapeModel.meanPoints(1:nPtsTibia,:));

    %Sort the registered points so they align with the reference mesh
    tibiaV_shapeModel = tibiaV(regSortIdxTib,:);

% % %     %Visualise to confirm that points now correspond. This can be done
% % %     %by using the target mesh faces against the registered moving
% % %     %points. The surface meshes should essentially come out looking the
% % %     %same (i.e. overlapped) despite using different faces
% % %     cFigure; hold on;
% % %     gpatch(tibiaF, tibiaV, 'rw', 'none', 0.3);
% % %     gpatch(tibTrabShapeModel.F1, tibiaV_shapeModel, 'gw', 'k');
% % %     title('Registered Tibia');
% % %     axisGeom;

    %Set a starting guess of the mean with appropriate bounds
    %Set the upper and lower bounds on the score variables. Here we set them as
    %-5 / +5 standard deviations about the PC score mean (i.e. zero)
    sdRange = 5;
    %Loop through retained PCs
    for nPC = 1:length(tibTrabShapeModel.varExplained)
        %Calculate scores and set to an initial guess variable
        x0(nPC,1) = 0;
        %Set the lower bound
        lb(nPC,1) = std(tibTrabShapeModel.score(:,nPC)) * -sdRange;
        %Set the upper bound
        ub(nPC,1) = std(tibTrabShapeModel.score(:,nPC)) * sdRange;
    end

    %Create the function handle for fmincon
    optFunc = @(pcScores)calcReconstructionError(pcScores, tibTrabShapeModel, tibiaV_shapeModel);
    
    %Set fmincon options
    fminconOpts = optimoptions('fmincon');
    fminconOpts.MaxIterations = 5000; %up potential iterations
    fminconOpts.MaxFunctionEvaluations = 20000; %up potential function evals
    fminconOpts.StepTolerance = 1e-4; %scale optimisation tolerance to problem
    fminconOpts.Display = 'iter'; %set to display function iterations

    %Run fmincon to optimise PC scores for trabecular reconstruction
    [x,fval] = fmincon(optFunc, x0, [], [], [], [], lb, ub, [], fminconOpts);

    %Reconstruct using the optimised PC scores
    optReconstructed = x(1:tibTrabShapeModel.retainPCs)' * ...
        tibTrabShapeModel.loadings(:,1:tibTrabShapeModel.retainPCs)' + ...
        tibTrabShapeModel.mean;

    %Reshape reconstructed points
    optReconstructedV = reshape(optReconstructed', [3, length(optReconstructed)/3])';
    
    %Extract the tibia and trabecular sections from reconstruction
    reconstructedV_tib = optReconstructedV(1:nPtsTibia,:);
    reconstructedV_trab = optReconstructedV(nPtsTibia+1:end,:);

% % %     %Visualise original vs. reconstructed
% % %     cFigure; hold on;
% % %     gpatch(tibTrabShapeModel.F1, tibiaV_shapeModel, 'gw', 'none', 0.3);
% % %     gpatch(tibTrabShapeModel.F1, reconstructedV_tib, 'rw', 'k', 1);
% % %     axisGeom; camlight headlight

% % %     %Visualise new predicted trabecular within reconstructed tibia
% % %     cFigure; hold on;
% % %     gpatch(tibTrabShapeModel.F1, reconstructedV_tib, 'gw', 'none', 0.3);
% % %     gpatch(tibTrabShapeModel.F2, reconstructedV_trab, 'bw', 'k', 1);
% % %     axisGeom; camlight headlight

    %Rigidly align the predicted trabecular to the original tibia surface
    [cpdRigTformTrabToOrig] = cpd_register(tibiaV_shapeModel, reconstructedV_trab, optCPD_rig);

    %Create new trabecular aligned points
    alignedNewTrabV = cpdRigTformTrabToOrig.Y;

% % %     %Visualise new predicted aligned trabecular within original tibia
% % %     cFigure; hold on;
% % %     gpatch(tibTrabShapeModel.F1, tibiaV_shapeModel, 'gw', 'none', 0.3);
% % %     gpatch(tibTrabShapeModel.F2, alignedNewTrabV, 'bw', 'k', 1);
% % %     axisGeom; camlight headlight

    %Identify trabecular points within (and hence also outside) of original tibia
    trabPtsIn = intriangulation(tibiaV_shapeModel, tibTrabShapeModel.F1, alignedNewTrabV, 0);
        
    %Rescale trabecular if not all points are inside tibia
    if sum(~trabPtsIn) > 0
        newPredictedTrabV = alignTrabecularModel(tibiaV_shapeModel, alignedNewTrabV, tibTrabShapeModel, trabPtsIn);
    end
    
% % %     %Produce final visualisation to check trabecular
% % %     cFigure; hold on;
% % %     gpatch(tibTrabShapeModel.F1, tibiaV_shapeModel, 'gw', 'none', 0.3);
% % %     gpatch(tibTrabShapeModel.F2, newPredictedTrabV, 'bw', 'k', 1);
% % %     axisGeom; camlight headlight
       
    %Allocate altered tibia faces and vertices to variables
    newTibiaF = tibTrabShapeModel.F1;
    newTibiaV = tibiaV_shapeModel;

end