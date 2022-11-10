function [trabV] = alignTrabecularModel(tibiaV, predictedV, shapeModel, trabPtsIn)

    %% This script aligns the generic trabecular model to the original tibia
    %
    %  Inputs:
    %       tibiaV - vertices for the tibia surface being used with prediction
    %       predictedV - vertices for the trabecular surface predicted
    %       shapeModel - the structure with shape model data for the tibia-plus-trabecular
    %       trabPtsIn - logical indicative of points inside and outside of tibia
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
        
    %Set options for remeshing
    optionStruct_tib.nb_pts = 3500; %Set desired number of points
    optionStruct_tib.disp_on = 0; % Turn off command window text display

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
    
    %Set default for smoothing
    if nargin < 5
        smoothIterations = 0;
    end

    %% Apply trabecular prediction to the input surface

    %Generate a generic trabecular to support an initial guess of rescaling 
    [genericTrabF, genericTrabV] = generateGenericTrabecular(...
        tibiaV, shapeModel, 1, length(predictedV));

    %Non-rigidly align the predicted trabecular to the generic trabecular
    %to identify matching points where corrections are needed
    [cpdRigTformTrab] = cpd_register(predictedV, genericTrabV, optCPD);

    %Identify the matching points in the target mesh against the
    %registration to identify the corresponding point indices
    regSortIdxTrab = knnsearch(cpdRigTformTrab.Y, predictedV);

    %Identify displacement of each point from predicted trabecular
    %Only extract displacement if the point in the tibia boundaries
    for ptNo = 1:length(predictedV)
        if ~trabPtsIn(ptNo)
            translationVals(ptNo,:) = predictedV(ptNo,:) - genericTrabV(regSortIdxTrab(ptNo),:);
        else
            translationVals(ptNo,:) = [0 0 0];
        end
    end

    %Create new trabecular vertices based on necessary translations
    alignedPredictedV = predictedV - translationVals;

    %One last final error check for if trabecular points are inside, as
    %later meshing functions will fail if this isn't the case
    trabPtsIn = intriangulation(tibiaV, shapeModel.F1, alignedPredictedV, 0);
    if sum(~trabPtsIn) > 0
        error('Something has gone wrong with containing trabecular points within tibia...whoops...')
    end
    
    %Check that any tibia points aren't outside the trabecular surface
    tibiaPtsIn = intriangulation(alignedPredictedV, shapeModel.F2, tibiaV, 0);
    
    %Correct if necessary
    if sum(tibiaPtsIn) > 0
        
        %Get point indices
        tibiaPtsInInd = find(tibiaPtsIn);
        
        %Loop through and correct points around
        for ptInd = 1:length(tibiaPtsInInd)
            
            %Get tibia point index as variable
            currTibiaPtInd = tibiaPtsInInd(ptInd);
            
            %Get the distance of the trabecular points to the current tibia point
            currPtDist = distancePoints3d(tibiaV(currTibiaPtInd,:), alignedPredictedV);
            
            %Find points less than the current 3mm radius
            manipulatePts = currPtDist <= 5;
            
            %Identify displacement of each point from predicted trabecular
            %Only extract displacement if the point in the tibia boundaries
            for ptNo = 1:length(alignedPredictedV)
                if manipulatePts(ptNo)
                    translationVals2(ptNo,:) = alignedPredictedV(ptNo,:) - cpdRigTformTrab.Y(regSortIdxTrab(ptNo),:);
                else
                    translationVals2(ptNo,:) = [0 0 0];
                end
            end
            
            %Replace predicted points by applying translations extracted
            alignedPredictedV = alignedPredictedV - translationVals2;
            
        end
  
    end
    
    %% Store data to output faces and vertices variables
    
    %Set outputs
    trabV = alignedPredictedV;

end