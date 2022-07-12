function [newPredictedTrabV, newTibiaF, newTibiaV] = createTrabecularFromSample(tibiaF, tibiaV, tibiaShapeModel, trabShapeModel, finalLinearModel)

    %% This script creates a trabecular surface from one of the sample
    %  tibia-fibula surfaces from the Nolte et al. dataset
    %
    %  Inputs:
    %       tibiaF - faces for the sample tibia surface
    %       tibiaV - vertices for the sample tibia surface
    %       tibiaShapeModel - the structure with shape model data for the tibia
    %       trabShapeModel - the structure with shape model data for the trabecular
    %       finalLinearModel - set of linear models for predicting PCs
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
    transMat = createTranslation3d(mean(tibiaShapeModel.meanPoints) - mean(tibiaV));
    %Apply translation to points
    tibiaV = transformPoint3d(tibiaV,transMat);

    %Apply a rigid transformation to match the new tibia points with the shape model
    %Don't apply the scaling from this rigid transform
    [cpdTformTibRig] = cpd_register(tibiaShapeModel.meanPoints, tibiaV, optCPD_rig);
    tibiaV = tibiaV * cpdTformTibRig.R' + repmat(cpdTformTibRig.t,1,length(tibiaV))';

    %To map the new surface to the shape model it's important to identify the
    %corresponding points on the surface to that of the shape model.

    %Perform non-rigid registration using CPD to align surfaces
    [cpdTformTib] = cpd_register(tibiaShapeModel.meanPoints, tibiaV, optCPD);

    %Identify the matching points in the target mesh against the
    %registration to identify the corresponding point indices
    regSortIdxTib = knnsearch(cpdTformTib.Y, tibiaShapeModel.meanPoints);

    %Sort the registered points so they align with the reference mesh
    tibiaV_shapeModel = tibiaV(regSortIdxTib,:);

    % % % %Visualise to confirm that points now correspond. This can be done
    % % % %by using the target mesh faces against the registered moving
    % % % %points. The surface meshes should essentially come out looking the
    % % % %same (i.e. overlapped) despite using different faces
    % % % cFigure; hold on;
    % % % gpatch(tibiaF, tibiaV, 'rw', 'none', 0.3);
    % % % gpatch(tibiaShapeModel.F, tibiaV_shapeModel, 'gw', 'k');
    % % % title('Registered Tibia');
    % % % axisGeom;

    %Reshape the registered data points and remove the mean
    newData = reshape(tibiaV_shapeModel', [1, length(tibiaV_shapeModel)*3]) - tibiaShapeModel.mean;

    %Project the new data against the loadings to find the estimated scores
    %Set the upper and lower bounds on the score variables. Here we set them as
    %-5 / +5 standard deviations about the PC score mean (i.e. zero)
    sdRange = 5;
    %Loop through retained PCs
    for nPC = 1:length(tibiaShapeModel.varExplained)
        %Calculate scores and set to an initial guess variable
        x0(nPC,1) = dot(tibiaShapeModel.loadings(:,nPC), newData);
        %Set the lower bound
        lb(nPC,1) = std(tibiaShapeModel.score(:,nPC)) * -sdRange;
        %Set the upper bound
        ub(nPC,1) = std(tibiaShapeModel.score(:,nPC)) * sdRange;
    end

    %Create the function handle for fmincon
    optFunc = @(pcScores)calcReconstructionError(pcScores, tibiaShapeModel, tibiaV_shapeModel);

    %Set fmincon options
    fminconOpts = optimoptions('fmincon');
    fminconOpts.MaxIterations = 5000; %up potential iterations
    fminconOpts.MaxFunctionEvaluations = 20000; %up potential function evals
    fminconOpts.StepTolerance = 1e-4; %scale optimisation tolerance to problem
    fminconOpts.Display = 'iter'; %set to display function iterations

    %Run fmincon to optimise PC scores for trabecular reconstruction
    [x,fval] = fmincon(optFunc, x0, [], [], [], [], lb, ub, [], fminconOpts);

    %Reconstruct using the optimised PC scores
    optReconstructed = x(1:tibiaShapeModel.retainPCs)' * ...
        tibiaShapeModel.loadings(:,1:tibiaShapeModel.retainPCs)' + ...
        tibiaShapeModel.mean;

    %Reshape reconstructed points
    optReconstructedV = reshape(optReconstructed', [3, length(optReconstructed)/3])';

    % % % %Visualise original vs. reconstructed
    % % % cFigure; hold on;
    % % % gpatch(tibiaShapeModel.F, tibiaV_shapeModel, 'gw', 'none', 0.3);
    % % % gpatch(tibiaShapeModel.F, optReconstructedV, 'rw', 'k', 1);
    % % % axisGeom; camlight headlight

    %Reconstruct a trabecular surface from previously created linear model
    for predictPC = 1:trabShapeModel.retainPCs

        %Output predicted values for left out case
        newTrabPCs(1, predictPC) = predict(finalLinearModel{predictPC}, x(1:tibiaShapeModel.retainPCs)');

    end

    %Reconstruct the new trabecular using the predicted PCs
    newTrabReconstructed = newTrabPCs(:,1:trabShapeModel.retainPCs) * ...
        trabShapeModel.loadings(:,1:trabShapeModel.retainPCs)';

    %Reshape the predicted points to 3D
    newTrabV = reshape(newTrabReconstructed(1,:) + trabShapeModel.mean, ...
            [3, length(trabShapeModel.mean)/3])';

    % % % %Visualise new predicted trabecular within reconstructed tibia
    % % % cFigure; hold on;
    % % % gpatch(tibiaShapeModel.F, optReconstructedV, 'gw', 'none', 0.3);
    % % % gpatch(trabShapeModel.F, newTrabV, 'bw', 'k', 1);
    % % % axisGeom; camlight headlight

    %Rigidly align the predicted trabecular to the original tibia surface
    [cpdRigTformTrabToOrig] = cpd_register(tibiaV_shapeModel, newTrabV, optCPD_rig);

    %Create new trabecular aligned points
    alignedNewTrabV = cpdRigTformTrabToOrig.Y;

    % % % %Visualise new predicted aligned trabecular within original tibia
    % % % cFigure; hold on;
    % % % gpatch(tibiaShapeModel.F, tibiaV_shapeModel, 'gw', 'none', 0.3);
    % % % gpatch(trabShapeModel.F, alignedNewTrabV, 'bw', 'k', 1);
    % % % axisGeom; camlight headlight

    %Identify trabecular points within (and hence also outside) of original tibia
    trabPtsIn = intriangulation(tibiaV_shapeModel, tibiaShapeModel.F, alignedNewTrabV, 0);

    %Rescale trabecular if not all points are inside tibia
    if sum(~trabPtsIn) > 0

        %Generate a generic trabecular to support an initial guess of rescaling 
        [genericTrabF, genericTrabV] = generateGenericTrabecular(...
            tibiaV_shapeModel, tibiaShapeModel, 1, length(alignedNewTrabV));

        %Non-rigidly align the predicted trabecular to the generic trabecular
        %to identify matching points where corrections are needed
        [cpdRigTformTrab] = cpd_register(genericTrabV, alignedNewTrabV, optCPD);

        %Identify the matching points in the target mesh against the
        %registration to identify the corresponding point indices
        regSortIdxTrab = knnsearch(cpdRigTformTrab.Y, alignedNewTrabV);

        %Identify displacement of each point from predicted trabecular
        %Only extract displacement if the point in the tibia boundaries
        for ptNo = 1:length(alignedNewTrabV)
            if ~trabPtsIn(ptNo)
                translationVals(ptNo,:) = alignedNewTrabV(ptNo,:) - cpdRigTformTrab.Y(regSortIdxTrab(ptNo),:);
            else
                translationVals(ptNo,:) = [0 0 0];
            end
        end

        %Create new trabecular vertices based on necessary translations
        alignedNewTrabV = alignedNewTrabV - translationVals;

        %Smooth out new trabecular surface if required
        %Note that this is the main reason how trabecular points can end up
        %outside the tibial surface. If this occurs then the will
        %progressively decrease smoothing iterations to zero until all
        %points lie within the boundaries
        smoothIterations = 2;
        while smoothIterations > 0

            %Set predicted points variable
            toSmoothV = alignedNewTrabV;

            %Run smoothing
            cPar.n = smoothIterations; %Number of iterations
            cPar.Method = 'LAP'; %Smooth method
            [newPredictedTrabV] = patchSmooth(trabShapeModel.F, toSmoothV, [], cPar);

            %Check points within tibia
            trabPtsIn = intriangulation(tibiaV, tibiaF, newPredictedTrabV, 0);
            if sum(~trabPtsIn) == 0
                %Break the while loop
                break
            else
                %Reduce smoothing iterations
                smoothIterations = smoothIterations - 1;
            end

        end

        %If smoothing iterations have progressed to zero then need to use
        %the original points
        if smoothIterations == 0
            newPredictedTrabV = alignedNewTrabV;    
        end

    else

        %Use original
        newPredictedTrabV = alignedNewTrabV;

    end
    
    %Check that any tibia points aren't outside the trabecular surface
    tibiaPtsIn = intriangulation(newPredictedTrabV, trabShapeModel.F, tibiaV, 0);
    
    %Correct if necessary
    if sum(tibiaPtsIn) > 0
        
        %Get point indices
        tibiaPtsInInd = find(tibiaPtsIn);
        
        %Loop through and correct points around
        for ptInd = 1:length(tibiaPtsInInd)
            
            %Get tibia point index as variable
            currTibiaPtInd = tibiaPtsInInd(ptInd);
            
            %Get the distance of the trabecular points to the current tibia point
            currPtDist = distancePoints3d(tibiaV(currTibiaPtInd,:), newPredictedTrabV);
            
            %Find points less than the current 5mm radius
            manipulatePts = currPtDist <= 5;
            
            %Identify displacement of each point from predicted trabecular
            %Only extract displacement if the point in the tibia boundaries
            for ptNo = 1:length(alignedNewTrabV)
                if manipulatePts(ptNo)
                    translationVals2(ptNo,:) = newPredictedTrabV(ptNo,:) - cpdRigTformTrab.Y(regSortIdxTrab(ptNo),:);
                else
                    translationVals2(ptNo,:) = [0 0 0];
                end
            end
            
            %Replace predicted points by applying translations extracted
            newPredictedTrabV = newPredictedTrabV - translationVals2;
            
        end
  
    end
    
    %Allocate altered tibia faces and vertices to variables
    newTibiaF = tibiaF;
    newTibiaV = tibiaV;

end