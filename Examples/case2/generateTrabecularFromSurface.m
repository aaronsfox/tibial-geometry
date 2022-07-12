%% This script provides a demonstration of how to generate a volumetric mesh
%  of the tibia trabecular by registering a different surface to the shape
%  model and then combining this with our regresion model for the
%  relationship between tibia surface shape and trabecular volume.
%
%  See the README.MD in this folder for more descriptive details on this
%  process and the source data.
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

%Set home directory
homeDir = pwd;

%Add supplementary code paths
addpath(genpath('supplementary'));

%Load the tibia and trabecular shape models required
load('..\..\ShapeModels\tibia\tibiaShapeModel.mat');
load('..\..\ShapeModels\trabecular-tibia\trabShapeModel.mat');

%Set array to rotate views of surfaces later in reconstruction process
surfaceRot = [-90, 0, 90, 180];

%Set list to label reconstruction subplots with
viewLabel = [{'Anterior'}, {'Lateral'}, {'Posterior'}, {'Medial'}];

%Navigate to segmentation directory to get case ID names
cd('..\..\Segmentation\');

%Grab the case names
f = dir();
for ff = 3:length(f)
    caseID{ff-2} = f(ff).name;
end

% % % %Set reference case number
% % % refCase = '147211';
% % % 
% % % %Reorder case list to make target mesh first
% % % caseID = [caseID(find(contains(caseID, refCase))), ...
% % %     caseID(1:find(contains(caseID, refCase))-1), ...
% % %     caseID(find(contains(caseID, refCase))+1:end)];

%Return to home directory
cd(homeDir);

%Set options for remeshing
optionStruct_tib.nb_pts = 3500; %Set desired number of points
optionStruct_tib.disp_on = 0; % Turn off command window text display
% % % optionStruct_fib.nb_pts = 2000; %Set desired number of points
% % % optionStruct_fib.disp_on = 0; % Turn off command window text display

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

%% Create a regression model between tibia and trabecular

%Here we create a series of regression models that use the retained tibia
%shape model scores to predict each of retained trabecular shape model
%scores. The theoretical basis for this is that the shape of the tibia can
%be effective in predicting the shape of the trabecular.
%
%The first step here is to ascertain the potential error by iterating
%through a leave-one-out cross validation approach. The series of linear
%models is fit while leaving one case out, and the values of that left out
%case predicted. We then reconstruct using the leave-out predictions and
%calculate the error between the original and predicted trabecular
%surfaces.

%Create an array that will sort through the leave one out approach
leaveOut = linspace(1, size(tibiaShapeModel.nodes,1), size(tibiaShapeModel.nodes,1));

%Loop through and predict trabecular PCs for left-out cases
for predictCase = 1:length(leaveOut)
    
    %Create a boolean mask to leave out the current case
    leaveOutBool = leaveOut ~= predictCase;
    
    %Get the retained PC scores from the tibia to use as predictors
    %Leave out appropriate case
    tibiaPCs = tibiaShapeModel.score(leaveOutBool,1:tibiaShapeModel.retainPCs);

    %Generate regression models for each of the retained trabecular PCs
    for predictPC = 1:trabShapeModel.retainPCs

        %Get the current scores as training for the regression model
        %Leave out appropriate case
        trabPCs = trabShapeModel.score(leaveOutBool,predictPC);
        
        %Fit the linear model on the data
        linearModel{predictPC} = fitlm(tibiaPCs, trabPCs);

        %Output predicted values for left out case
        predictTrabPCs(predictCase, predictPC) = predict(linearModel{predictPC}, ...
            tibiaShapeModel.score(predictCase,1:tibiaShapeModel.retainPCs));

    end
    
end

%Create a scatter subplot of the predicted trabecular PCs
f = figure; f.Position = [100 100 1080 800];
for predictPC = 1:trabShapeModel.retainPCs
    
    %Create subplot
    subplot(2,2,predictPC); hold on
    
    %Plot data
    scatter(trabShapeModel.score(:,predictPC), predictTrabPCs(:,predictPC), ...
        'k', 'filled');
    
    %Plot line of equality
    ax = gca;
    plot([ax.XLim(1), ax.XLim(2)], [ax.XLim(1), ax.XLim(2)], 'r--');
    
    %Add title
    title(['Trabecular Shape Model PC',num2str(predictPC)]);
    
    %Add labels
    xlabel('Actual Shape Model Score');
    ylabel('Predicted Shape Model Score');
    
end

%%%%% TODO: save this scatter with details somewhere...?

%Reconstruct each case using the predicted PCs
predictReconstructed = predictTrabPCs(:,1:trabShapeModel.retainPCs) * ...
    trabShapeModel.loadings(:,1:trabShapeModel.retainPCs)';

%Set-up waitbar for processing trabecular surfaces
wbar_trab = waitbar(0, 'Predicting trabecular surfaces...');

%Loop through predicted reconstructions, calculate the mean and peak error,
%and visualise the reconstruction accuracy
for predictCase = 1:10%%%%% FIX ME!!! length(leaveOut)
    
    %Get the original surface in XYZ format
    actualV = reshape(trabShapeModel.nodes(predictCase,:), ...
        [3, length(trabShapeModel.mean)/3])';
    
    %Get the predicted reconstruction in XYZ format
    predictedV = reshape(predictReconstructed(predictCase,:) + trabShapeModel.mean, ...
        [3, length(trabShapeModel.mean)/3])';
    
    %Get the tibial surface for alignment
    tibiaV = reshape(tibiaShapeModel.nodes(predictCase,:), ...
        [3, length(tibiaShapeModel.mean)/3])';
       
    %Check need to scale trabecular to fit inside tibia
    
    %Determine trabecular points inside vs. outside tibia
    trabPtsIn = intriangulation(tibiaV, tibiaShapeModel.F, predictedV, 0);
    
    %Extract the initial number of points outside        
    trabPointsOutsideInitial(predictCase,1) = sum(~trabPtsIn);
    
    %Rescale trabecular if not all points are inside tibia
    if sum(~trabPtsIn) > 0
        
        %Generate a generic trabecular to support an initial guess of rescaling 
        [genericTrabF, genericTrabV] = generateGenericTrabecular(...
            tibiaV, tibiaShapeModel, 1, length(predictedV));
        
        %Non-rigidly align the predicted trabecular to the generic trabecular
        %to identify matching points where corrections are needed
% % %         [cpdRigTformTrab] = cpd_register(genericTrabV, predictedV, optCPD_rig);
        [cpdRigTformTrab] = cpd_register(genericTrabV, predictedV, optCPD);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Non-rigid approach

        %Identify the matching points in the target mesh against the
        %registration to identify the corresponding point indices
        regSortIdxTrab = knnsearch(cpdRigTformTrab.Y, predictedV);
        
        %Identify initial guess displacement of each point from predicted trabecular
        %Only extract an actual initial movement guess if the point wasn't
        %originally in the tibia boundaries
        for ptNo = 1:length(predictedV)
            if ~trabPtsIn(ptNo)
                translationVals(ptNo,:) = predictedV(ptNo,:) - cpdRigTformTrab.Y(regSortIdxTrab(ptNo),:);
            else
                translationVals(ptNo,:) = [0 0 0];
            end
        end
        
        %Create new trabecular vertices based on necessary translations
        newV = predictedV - translationVals;
        
        %Smooth out new trabecular surface
        cPar.n = 2; %Number of iterations
        cPar.Method = 'LAP'; %Smooth method
        [newV] = patchSmooth(trabShapeModel.F, newV, [], cPar);
        
        
        %%%%%%%%% TODO: clean up and use this newV value as the predictedV
        %%%%%%%%% value for subsequent calculations - seems to be the best
        %%%%%%%%% way to end up with the relevant shape and maintaining
        %%%%%%%%% points --- will need to test on subsequent non shape
        %%%%%%%%% model tibias later in script...
        
        
        
% % %         %Reshape to a singular vector
% % %         x0 = [translationVals(:,1); translationVals(:,2); translationVals(:,3)];
% % %         
% % %         %Create the function handle for fmincon
% % %         optFuncTrab = @(x)calcTrabecularError_v4(x, predictedV);
% % % 
% % %         %Create the constraint function handle for fmincon
% % %         conFuncTrab = @(x)nTrabPointsOutside_v4(x, predictedV, tibiaV, tibiaShapeModel);
% % %       
% % %         %Set upper and lower bounds for point movement
% % %         xLB = ones(length(x0),1) * -10;
% % %         xUB = ones(length(x0),1) * 10;
% % %         
% % %         %Set fmincon options
% % %         fminconOpts = optimoptions('fmincon');
% % %         fminconOpts.MaxIterations = 5000; %up potential iterations
% % %         fminconOpts.MaxFunctionEvaluations = 20000; %up potential function evals
% % %         fminconOpts.StepTolerance = 1e-6; %scale optimisation tolerance to problem
% % %         fminconOpts.Display = 'iter'; %set to display function iterations
% % %         
% % %         %Run fmincon
% % %         [x_trab, fval_trab] = fmincon(optFuncTrab, x0, [], [], [], [], xLB, xUB, conFuncTrab, fminconOpts);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%

        
        

% % %         %Sort the registered points so they align with the reference mesh
% % %         tibiaV_shapeModel = tibiaV(regSortIdxTib,:);
        
% % %         %Extract relevant factors as initial guess for optimisation
% % %         x0 = reshape([cpdRigTformTrab.s, cpdRigTformTrab.s, cpdRigTformTrab.s; ... %apply generic scaling over axes
% % %             rotm2eul(cpdRigTformTrab.R, 'XYZ'); ... %convert rotation matrix to XYZ Euler sequence
% % %             cpdRigTformTrab.t']', ... %extract translation factors
% % %             [1,9]);
% % %         
% % %         %Set the weights for different parameters in the optimisation
% % %         w = [0.1, 10, 0.1, ... %weighting for XYZ scaling
% % %             2.5, 0.1, 2.5, ... %weighting for XYZ rotation
% % %             0.01, 0.1, 0.01, ... %weighting for XYZ translation
% % %             100]; %weighting for number of points outside
% % %         
% % %         %Create the function handle for fmincon
% % %         optFuncTrab = @(x)optimiseTrabecularPointsOutside(x, predictedV, tibiaV, tibiaShapeModel, w);
% % % 
% % % % % %         %Create the constraint function handle for fmincon
% % % % % %         conFuncTrab = @(x)nTrabPointsOutside_v3(x, predictedV, tibiaV, tibiaShapeModel);  
% % % 
% % %         %Set upper and lower bounds for scale factors
% % %         %Lower bounds
% % %         xLB = [0.7, 0.7, 0.7, ... %XYZ scaling
% % %             -3, -10, -3, ... %XYZ rotation
% % %             -5, -20, -5]; %XYZ translation
% % %         %Upper bounds
% % %         xUB = [1, 1, 1, ... %XYZ scaling
% % %             3, 10, 3, ... %XYZ rotation
% % %             5, 20, 5]; %XYZ translation
% % % 
% % %         %Set fmincon options
% % %         fminconOpts = optimoptions('fmincon');
% % %         fminconOpts.MaxIterations = 5000; %up potential iterations
% % %         fminconOpts.MaxFunctionEvaluations = 20000; %up potential function evals
% % %         fminconOpts.StepTolerance = 1e-6; %scale optimisation tolerance to problem
% % %         fminconOpts.Display = 'iter'; %set to display function iterations
% % %         
% % %         %Run fmincon
% % %         [x_trab] = fmincon(optFuncTrab, x0, [], [], [], [], xLB, xUB, [], fminconOpts);
% % % 
% % %         %Reconstruct with optimised scale factors
% % %         scaledPredictedV = x_trab(1:3) .* predictedV;
% % %         
% % %         %Recentre trabecular to tibia   
% % %         %Get central points along both point clouds
% % %         centralTibia = mean([max(tibiaV);min(tibiaV)]);
% % %         centralTrab = mean([max(scaledPredictedV);min(scaledPredictedV)]);    
% % % 
% % %         %Adjust the predicted trabecular to centre of tibia
% % %         scaledPredictedV = scaledPredictedV - (centralTrab - centralTibia);
% % %         
% % %         %Apply the optimised rotational values
% % %         scaledPredictedV = transformPoint3d(scaledPredictedV, ...
% % %             createRotationOx(deg2rad(x_trab(4))));
% % %         scaledPredictedV = transformPoint3d(scaledPredictedV, ...
% % %             createRotationOy(deg2rad(x_trab(5))));
% % %         scaledPredictedV = transformPoint3d(scaledPredictedV, ...
% % %             createRotationOz(deg2rad(x_trab(6))));
% % %         
% % %         %Apply the optimised translation values
% % %         scaledPredictedV = scaledPredictedV - x_trab(7:9);
% % %         
        %Check number of points outside
        scaledPtsIn = intriangulation(tibiaV, tibiaShapeModel.F, newV, 0);
% % %         
% % %         %Extract the final number of points outside        
% % %         trabPointsOutsideFinal(predictCase,1) = sum(~scaledPtsIn);
% % %         
        %Reset predicted trabecular points
        predictedV = newV;
                
    else
        
        %Set trabecular points outside variable to zero
        trabPointsOutsideFinal(predictCase,1) = 0;
        
    end
    
    %Calculate the Jaccard index between actual and predicted trabecular
    jaccardSimilarity(predictCase,1) = calcJaccardTrabecular(actualV, predictedV, trabShapeModel);
    
    %Calculate point error distance
    pointErrorDist(predictCase,:) = distancePoints3d(actualV, predictedV);

    %Calculate the mean error
    pointErrorDistMean(predictCase,1) = mean(pointErrorDist(predictCase,:));
    
    %Calculate the peak error
    pointErrorDistMax(predictCase,1) = max(pointErrorDist(predictCase,:));
    
    %Convert distance error to colour scales for visualisation
    pointErrorDistColF(predictCase,:) = vertexToFaceMeasure(trabShapeModel.F, pointErrorDist(predictCase,:)');
    pointErrorDistColV(predictCase,:) = faceToVertexMeasure(trabShapeModel.F, actualV, pointErrorDistColF(predictCase,:)');
    
% % %     %Create visualisation of error
% % %     %Use subplots to create different perspectives
% % %     cFigure; hold on;
% % %     %Loop through four views to create subplot
% % %     for viewNo = 1:4    
% % %         %Create subplot for current view
% % %         subplot(1,4,viewNo);
% % %         %Add greyed out original surface
% % %         hpOrig = gpatch(trabShapeModel.F, actualV, [200/255 200/255 200/255], 'none', 0.5);
% % %         %Add colormapped reconstruction
% % %         hpPred = gpatch(trabShapeModel.F, predictedV, pointErrorDistColV(predictCase,:)', 'none', 1);
% % %         %Interpolate colouring for smoothness
% % %         hpPred.FaceColor = 'Interp'; colormap viridis
% % %         %Set axis view
% % %         axis equal; axis tight; view(0,90);
% % %         rotate(hpOrig,[0 1 0], surfaceRot(viewNo));
% % %         rotate(hpPred,[0 1 0], surfaceRot(viewNo));
% % %         %Set axis parameters
% % %         camlight headlight; axis off
% % %         %Add colorbar on last view
% % %         if viewNo == 4
% % %             colorbar
% % %         end
% % %         %Add title
% % %         title([viewLabel{viewNo},' View'], 'FontSize', 12);
% % %     end
% % %     
% % %     %Export figure
% % %     export_fig(['figures\predictedReconstructionError\',caseID{predictCase},'_predictedReconstructionErrorMap.png'],'-m1');
% % % 	close

    %Update waitbar
    waitbar(predictCase/length(leaveOut), wbar_trab);
   
end

%Calculate and display mean and 95% CI's for the mean and peak error
%Mean & SD for mean error
meanError_m = mean(pointErrorDistMean);
meanError_sd = std(pointErrorDistMean);
meanError_lower95 = meanError_m - (1.96 * (meanError_sd / sqrt(length(leaveOut))));
meanError_upper95 = meanError_m + (1.96 * (meanError_sd / sqrt(length(leaveOut))));
%Mean & SD for maxerror
maxError_m = mean(pointErrorDistMax);
maxError_sd = std(pointErrorDistMax);
maxError_lower95 = maxError_m - (1.96 * (maxError_sd / sqrt(length(leaveOut))));
maxError_upper95 = maxError_m + (1.96 * (maxError_sd / sqrt(length(leaveOut))));
%Display
disp(['Reconstruction mean point error (mean +/- 95% CIs) = ',num2str(round(meanError_m,2)), ...
    '[',num2str(round(meanError_lower95,2)),',',num2str(round(meanError_upper95,2)),']']);
disp(['Reconstruction max point error (mean +/- 95% CIs) = ',num2str(round(maxError_m,2)), ...
    '[',num2str(round(maxError_lower95,2)),',',num2str(round(maxError_upper95,2)),']']);

%Create the final linear models using all cases

%Get the retained PC scores from the tibia to use as predictors
%Leave out appropriate case
tibiaPCs = tibiaShapeModel.score(:,1:tibiaShapeModel.retainPCs);

%Generate regression models for each of the retained trabecular PCs
for predictPC = 1:trabShapeModel.retainPCs

    %Get the current scores as training for the regression model
    %Leave out appropriate case
    trabPCs = trabShapeModel.score(:,predictPC);

    %Fit the linear model on the data
    finalLinearModel{predictPC} = fitlm(tibiaPCs, trabPCs);
    
% % %     %Display model output
% % %     finalLinearModel{predictPC}

% % %     %Plot model
% % %     figure; plot(finalLinearModel{predictPC});

end

%% Apply trabecular prediction to a new surface

%%%%% Scaling down trabecular to be within tibia is troubling --- idea with
%%%%% equality constraints in an optimisation seems right...but just isn't
%%%%% working as expected...think about it a little more...

%%%%% TODO: consider using another of the NMDID database surfaces so that
%%%%% the population matches appropriately...

%Load the tibia and fibula surface meshes
% % % %Tibia
% % % [tibiaSTLstruct] = import_STL('sampleSurfaces\rightTibia.stl');
% % % tibiaF = tibiaSTLstruct.solidFaces{1}; %Faces
% % % tibiaV = tibiaSTLstruct.solidVertices{1}; %Vertices
% % % [tibiaF,tibiaV] = mergeVertices(tibiaF,tibiaV);
% % % %Fibula
% % % [fibulaSTLstruct] = import_STL('sampleSurfaces\rightFibula.stl');
% % % fibulaF = fibulaSTLstruct.solidFaces{1}; %Faces
% % % fibulaV = fibulaSTLstruct.solidVertices{1}; %Vertices
% % % [fibulaF,fibulaV] = mergeVertices(fibulaF,fibulaV);

%Tibia-fibula
[tibFibSTLstruct] = import_STL('sampleSurfaces\Nolte2016\C06RTF.stl');
tibFibF = tibFibSTLstruct.solidFaces{1}; %Faces
tibFibV = tibFibSTLstruct.solidVertices{1}; %Vertices
[tibFibF,tibFibV] = mergeVertices(tibFibF,tibFibV);

%Split the surfaces to grab the tibia only

%Group the vertices and faces from the combined shape model
[groupIndexVertices,groupIndexFaces] = groupVertices(tibFibF, tibFibV,0);

% % % %Visualise
% % % cFigure; hold on
% % % title('Grouped Faces')
% % % gpatch(tibFibF, tibFibV, groupIndexFaces,'none');
% % % axisGeom; camlight headlight;
% % % colormap gjet; icolorbar;

%Identify which grouped section contains a higher volume
%This will be indicative of the tibia
if tetVolMeanEst(tibFibF(groupIndexFaces == 1,:),tibFibV) > ...
        tetVolMeanEst(tibFibF(groupIndexFaces == 2,:),tibFibV)
    %First index is tibia
    logicKeep = groupIndexFaces == 1;
else
    %Second index is tibia
    logicKeep = groupIndexFaces == 2;
end
%Separate the surfaces
[tibiaF, tibiaV] = patchCleanUnused(tibFibF(logicKeep,:), tibFibV);
[fibulaF, fibulaV] = patchCleanUnused(tibFibF(~logicKeep,:), tibFibV);

% % % %Visualise
% % % cFigure; hold on
% % % gpatch(tibiaF,tibiaV,'gw', 'none')
% % % gpatch(fibulaF,fibulaV,'bw', 'none')
% % % axisGeom; camlight headlight

%Remesh the surfaces to match the shape model points
[tibiaF,tibiaV] = ggremesh(tibiaF, tibiaV, optionStruct_tib);

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

%Apply a simple procrustes transform to maximise chance of best alignment
[~, tibiaV] = procrustes(tibiaShapeModel.meanPoints, tibiaV, 'scaling', false);

%To map the new surface to the shape model it's important to identify the
%corresponding points on the surface to that of the shape model.

%Perform non-rigid registration using CPD to align surfaces
[cpdTformTib] = cpd_register(tibiaShapeModel.meanPoints, tibiaV, optCPD);

%Identify the matching points in the target mesh against the
%registration to identify the corresponding point indices
regSortIdxTib = knnsearch(cpdTformTib.Y, tibiaShapeModel.meanPoints);

%Sort the registered points so they align with the reference mesh
tibiaV_shapeModel = tibiaV(regSortIdxTib,:);

%Visualise to confirm that points now correspond. This can be done
%by using the target mesh faces against the registered moving
%points. The surface meshes should essentially come out looking the
%same (i.e. overlapped) despite using different faces
% % % cFigure; hold on;
% % % gpatch(tibiaF, tibiaV, 'rw', 'none', 0.3);
% % % gpatch(tibiaShapeModel.F, tibiaV_shapeModel, 'gw', 'k');
% % % title('Registered Tibia');
% % % axisGeom;

% % % %Apply procrustes transform to rigidly align surfaces
% % % [~, regTibV, procTransform] = procrustes(tibiaShapeModel.meanPoints, ...
% % %     tibiaV_shapeModel, 'scaling', false);

% % % %Apply rigid transform to rigidly align surfaces
% % % [cpdRigTformTib] = cpd_register(tibiaShapeModel.meanPoints, tibiaV_shapeModel, optCPD_rig);

% % % %Apply rigid transform to points
% % % %Combine rotation/transformation matrix
% % % rigTransform = [cpdRigTformTib.R, cpdRigTformTib.t];
% % % % % % %Combine rotation/transformation/scaling matrix
% % % % % % rigTransform = [cpdRigTformTib.R, cpdRigTformTib.t; 0, 0, 0, cpdRigTformTib.s];
% % % %Transform points
% % % regTibV(ptNo,:) = transformPoint3d(tibiaV_shapeModel(ptNo,:), rigTransform);
%Direct application of rotation, transformation and scaling
% % % regTibV = cpdRigTformTib.s * tibiaV_shapeModel * cpdRigTformTib.R' + ...
% % %     repmat(cpdRigTformTib.t',[length(tibiaV_shapeModel) 1]);

%%%% TODO: this rigid transform seems to align to the base of the tibia
%%%% when the scaling isn't applied...

% % % %Visualise effect of rigid transform
% % % cFigure; hold on
% % % plotV(tibiaShapeModel.meanPoints, 'go');
% % % plotV(tibiaV_shapeModel, 'r.');
% % % plotV(regTibV, 'b.');
% % % % % % plotV(cpdRigTformTib.Y, 'm.');
% % % axisGeom;
% % % legend('Shape Model Mean Target', 'New Surface Points', 'New Surface Transformed Points');

% % % %Extract rigidly registered points
% % % regTibV = cpdRigTformTib.Y;

%%%% TODO: the score projection code below doesn't really work --- need to
%%%% test out fmincon approach that minimises the mean point error using
%%%% the first 5 PCs. An addition to this might be the need to optimise a
%%%% rigid transformation matrix that is also exported, so that when we
%%%% calculate the trabecular from the estimated PC scores we can then
%%%% apply that rigid transformation so that the trabecular aligns well
%%%% with the original bone surface...

%Generate an appropriate starting guess for the function by projecting the
%new surface data point on to the loadings of the shape model

%Reshape the registered data points and remove the mean
newData = reshape(tibiaV_shapeModel', [1, length(tibiaV_shapeModel)*3]) - tibiaShapeModel.mean;
% % % newData = reshape(regTibV', [1, length(regTibV)*3]) - tibiaShapeModel.mean;

%Project the new data against the loadings to find the estimated scores
%Set the upper and lower bounds on the score variables. Here we set them as
%-5 / +5 standard deviations about the PC score mean (i.e. zero)
sdRange = 5;
%Loop through retained PCs
% % % for nPC = 1:tibiaShapeModel.retainPCs
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
% % % optFunc = @(pcScores)calcReconstructionError(pcScores, tibiaShapeModel, regTibV);

%Set fmincon options
fminconOpts = optimoptions('fmincon');
fminconOpts.MaxIterations = 5000; %up potential iterations
fminconOpts.MaxFunctionEvaluations = 20000; %up potential function evals
fminconOpts.StepTolerance = 1e-4; %scale optimisation tolerance to problem
fminconOpts.Display = 'iter'; %set to display function iterations

%Run fmincon
[x,fval] = fmincon(optFunc, x0, [], [], [], [], lb, ub, [], fminconOpts);
% % % [x,fval] = fmincon(optFunc, zeros(5,1), [], [], [], [], lb, ub);
% % % calcReconstructionError(pcScores, tibiaShapeModel, regTibV);
% % % meanError = calcReconstructionError(x0, tibiaShapeModel, regTibV)

%Reconstruct using initial guess PC scores
% % % igReconstructed = x0(1:tibiaShapeModel.retainPCs)' * ...
% % %     tibiaShapeModel.loadings(:,1:tibiaShapeModel.retainPCs)' + ...
% % %     tibiaShapeModel.mean;
% % % igReconstructed = x0(1:length(tibiaShapeModel.varExplained))' * ...
% % %     tibiaShapeModel.loadings(:,1:length(tibiaShapeModel.varExplained))' + ...
% % %     tibiaShapeModel.mean;

%Reconstruct using the optimised PC scores
optReconstructed = x(1:tibiaShapeModel.retainPCs)' * ...
    tibiaShapeModel.loadings(:,1:tibiaShapeModel.retainPCs)' + ...
    tibiaShapeModel.mean;

%Reshape reconstructed points
% % % igReconstructedV = reshape(igReconstructed', [3, length(igReconstructed)/3])';
optReconstructedV = reshape(optReconstructed', [3, length(optReconstructed)/3])';

% % % %Visualise original vs. reconstructed
% % % cFigure; hold on;
% % % gpatch(tibiaShapeModel.F, tibiaV_shapeModel, 'gw', 'none', 0.3);
% % % % % % gpatch(tibiaShapeModel.F, regTibV, 'gw', 'k', 0.3);
% % % gpatch(tibiaShapeModel.F, optReconstructedV, 'rw', 'k', 1);
% % % % % % gpatch(tibiaShapeModel.F, igReconstructedV, 'bw', 'k', 1);
% % % axisGeom; camlight headlight

% % % %Apply procrustes algorithm between reconstruction and original
% % % [D, Z, TRANSFORM] = procrustes(tibiaV_shapeModel, optReconstructedV, 'scaling', true);

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
% % % % % % gpatch(tibiaShapeModel.F, tibiaV_shapeModel, 'rw', 'none', 0.3);
% % % % % % gpatch(tibiaShapeModel.F, regTibV, 'gw', 'none', 0.3);
% % % gpatch(tibiaShapeModel.F, optReconstructedV, 'gw', 'none', 0.3);
% % % gpatch(trabShapeModel.F, newTrabV, 'bw', 'k', 1);
% % % axisGeom; camlight headlight

%Rigidly align the predicted trabecular to the original tibia surface
[cpdRigTformTrabToOrig] = cpd_register(tibiaV_shapeModel, newTrabV, optCPD_rig);

%Create new trabecular aligned points
alignedNewTrabV = cpdRigTformTrabToOrig.Y;
% % % 
% % % %Visualise new predicted and aligned trabecular within original tibia
% % % cFigure; hold on;
% % % gpatch(tibiaShapeModel.F, tibiaV_shapeModel, 'gw', 'none', 0.3);
% % % gpatch(trabShapeModel.F, alignedNewTrabV, 'bw', 'k', 1);
% % % axisGeom; camlight headlight

%Identify trabecular points within (and hence also outside) of tibia
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

    %Smooth out new trabecular surface
    cPar.n = 2; %Number of iterations
    cPar.Method = 'LAP'; %Smooth method
    [newV] = patchSmooth(trabShapeModel.F, alignedNewTrabV, [], cPar);

else
    
    %Use original
    newV = alignedNewTrabV;
    
end

%Review newly fixed trabecular within original tibia
cFigure; hold on;
gpatch(tibiaShapeModel.F, tibiaV_shapeModel, 'gw', 'none', 0.3);
% % % gpatch(genericTrabF, genericTrabV, 'bw', 'k', 1);
gpatch(trabShapeModel.F, newV, 'bw', 'k', 1);
axisGeom; camlight headlight

% % % %Check for new points in
% % % %%%% Still some potential for minimal number of points to sit outside...
% % % trabPtsIn = intriangulation(tibiaV_shapeModel, tibiaShapeModel.F, newV, 0);
% % % sum(~trabPtsIn);

%%%%%%%% BELOW IS OLD TESTING...

% % % %Set a scaled copy of the reconstructed trabecular if necessary
% % % if ~all(trabPtsIn)
% % %     newTrabV_scaled = newTrabV;
% % % end

%%%%% This approach doesn't really work as it just makes it tiny
% % % %Check if there are points outside
% % % while ~all(trabPtsIn)
% % %     
% % %     %Scale down trabecular by 1% increment
% % %     newTrabV_scaled = newTrabV_scaled * 0.99;
% % %     
% % %     %Re-allign using procrustes analysis
% % %     [~, newTrabV_scaled] = procrustes(newTrabV, newTrabV_scaled, 'scaling', false);
% % %     
% % %     %Re-check if all trab points are inside
% % %     trabPtsIn = intriangulation(optReconstructedV, tibiaShapeModel.F, newTrabV_scaled, 0);
% % %     
% % % end
%%%%% Need to slightly morphy trabecular PCs until surface is inside the
%%%%% tibia


%%%%% This optimisation process seems like it should work if the parameters
%%%%% get right...Something about the bounds that is making the
%%%%% optimisation act in a weird way, or the initial guess skewing the
%%%%% final result...


% % % %Create the function handle for fmincon
% % % optFuncTrab = @(x)calcTrabecularError(x, trabShapeModel, newTrabV);
% % % % % % optFuncTrab = @(x)calcTrabecularError_v2(x, trabShapeModel, newTrabV);
% % % 
% % % %Create the constraint function handle for fmincon
% % % conFuncTrab = @(x)nTrabPointsOutside(x, trabShapeModel, optReconstructedV, tibiaShapeModel.F);
% % % % % % conFuncTrab = @(x)nTrabPointsOutside_v2(x, newTrabV, optReconstructedV, tibiaShapeModel.F);
% % % 
% % % % % % %Starting scaling factor
% % % % % % scaleFactor = [1, 1, 1];
% % % % % % 
% % % % % % %Starting translation matrix
% % % % % % transMat = [0, 0, 0];
% % % % % % 
% % % % % % %Starting rotation matrix
% % % % % % rotMat = [1, 0, 0, 0, 1, 0, 0, 0, 1];
% % % % % % 
% % % % % % x0_trab = [newTrabPCs, scaleFactor, transMat, rotMat];
% % % % % % x0_trab = [0, 0, 0, 0, scaleFactor, transMat, rotMat];
% % % 
% % % % % % %Set the starting guess for the optimisation function as the predicted
% % % % % % %trabecular PC scores
% % % % % % x0_trab = newTrabPCs;
% % % % % % 
% % % % % % %Set the boundaries to search around the predicted trabecular using
% % % % % % sdRange_trab = 3;
% % % % % % %Loop through retained PCs
% % % % % % for nPC = 1:length(x0_trab)
% % % % % %     %Set the lower bound
% % % % % %     lb_trab(1,nPC) = x0_trab(nPC) - (std(trabShapeModel.score(:,nPC)) * sdRange_trab);
% % % % % %     %Set the upper bound
% % % % % %     ub_trab(1,nPC) = x0_trab(nPC) + (std(trabShapeModel.score(:,nPC)) * sdRange_trab);
% % % % % % end
% % % 
% % % sdRange = 3;
% % % %Loop through retained PCs
% % % % % % for nPC = 1:tibiaShapeModel.retainPCs
% % % for nPC = 1:length(newTrabPCs)
% % %     %Calculate scores and set to an initial guess variable
% % %     x0_trab(1,nPC) = 0;
% % %     %Set the lower bound
% % %     lb_trab(1,nPC) = std(trabShapeModel.score(:,nPC)) * -sdRange;
% % %     %Set the upper bound
% % %     ub_trab(1,nPC) = std(trabShapeModel.score(:,nPC)) * sdRange;
% % % end
% % % 
% % % %Other settings and bounds
% % % %Scale factors
% % % % % % x0_trab = [x0_trab, 1, 1, 1];
% % % x0_trab = [x0_trab, 0.97, 0.95, 0.97];
% % % lb_trab = [lb_trab, 0.9, 0.8, 0.9];
% % % ub_trab = [ub_trab, 1.1, 1.2, 1.1];
% % % %Translation matrix
% % % x0_trab = [x0_trab, 0, 0, 0];
% % % lb_trab = [lb_trab, -50, -300, -20];
% % % ub_trab = [ub_trab, 50, 300, 50];
% % % %Rotation matrix
% % % x0_trab = [x0_trab, 1, 0, 0, 0, 1, 0, 0, 0, 1];
% % % lb_trab = [lb_trab, -pi, -pi, -pi, -pi, -pi, -pi, -pi, -pi, -pi];
% % % ub_trab = [ub_trab, pi, pi, pi, pi, pi, pi, pi, pi, pi];
% % % 
% % % %Set fmincon options
% % % fminconOpts = optimoptions('fmincon');
% % % fminconOpts.MaxIterations = 5000; %up potential iterations
% % % fminconOpts.MaxFunctionEvaluations = 20000; %up potential function evals
% % % fminconOpts.StepTolerance = 1e-8; %scale optimisation tolerance to problem
% % % 
% % % %Run fmincon
% % % [x_trab,fval_trab] = fmincon(optFuncTrab, x0_trab, [], [], [], [], lb_trab, ub_trab, conFuncTrab, fminconOpts);
% % % 
% % % %Extract the shape model scores from the function output
% % % optTrabPCs(1:trabShapeModel.retainPCs,1) = x_trab(1:trabShapeModel.retainPCs)';
% % % 
% % % %Extract the scaling factor
% % % optScaleFactor(1:3,1) = x_trab(trabShapeModel.retainPCs+1:trabShapeModel.retainPCs+3);
% % % % % % optScaleFactor(1:3,1) = x_trab(1:3);
% % % 
% % % %Extract the translational matrix
% % % optTransMat(1:3,1) = x_trab(trabShapeModel.retainPCs+4:trabShapeModel.retainPCs+6);
% % % % % % optTransMat(1:3,1) = x_trab(4:6);
% % % 
% % % %Extract the rotation matrix
% % % optRotMat(1,1:3) = x_trab(trabShapeModel.retainPCs+7:trabShapeModel.retainPCs+9);
% % % optRotMat(2,1:3) = x_trab(trabShapeModel.retainPCs+10:trabShapeModel.retainPCs+12);
% % % optRotMat(3,1:3) = x_trab(trabShapeModel.retainPCs+13:trabShapeModel.retainPCs+15);
% % % % % % optRotMat(1,1:3) = x_trab(7:9);
% % % % % % optRotMat(2,1:3) = x_trab(10:12);
% % % % % % optRotMat(3,1:3) = x_trab(13:15);
% % % 
% % % %Reconstruct the new trabecular using the optimised PCs
% % % newTrabOptimised = optTrabPCs' * ...
% % %     trabShapeModel.loadings(:,1:trabShapeModel.retainPCs)';
% % % % % % newTrabOptimised = newTrabPCs * ...
% % % % % %     trabShapeModel.loadings(:,1:trabShapeModel.retainPCs)';
% % %    
% % % %Reshape the predicted points to 3D
% % % optTrabV = reshape(newTrabOptimised(1,:) + trabShapeModel.mean, ...
% % %         [3, length(trabShapeModel.mean)/3])';
% % %     
% % % %Apply translation, rotation and scaling
% % % finalTrabV = optScaleFactor' .* optTrabV * optRotMat' + ...
% % %     repmat(optTransMat',[length(optTrabV) 1]);
% % % 
% % % %Identify trabecular points within (and hence also outside) of tibia
% % % trabPtsIn = intriangulation(optReconstructedV, tibiaShapeModel.F, finalTrabV, 0);
% % % 
% % % %%%% UP TO WITHIN HERE FOR TESTING... %%%%%
% % % 
% % % %Visualise new predicted trabecular w/ outside points
% % % cFigure; hold on;
% % % % % % gpatch(tibiaShapeModel.F, tibiaV_shapeModel, 'gw', 'none', 0.3);
% % % % % % gpatch(tibiaShapeModel.F, regTibV, 'gw', 'none', 0.3);
% % % gpatch(tibiaShapeModel.F, optReconstructedV, 'gw', 'none', 0.3);
% % % % % % gpatch(trabShapeModel.F, newTrabV, 'bw', 'k', 1);
% % % gpatch(trabShapeModel.F, finalTrabV, 'rw', 'k', 1);
% % % plotV(finalTrabV(~trabPtsIn,:), 'b.');
% % % axisGeom; camlight headlight
% % % 
% % % %%%%% NOTE: predicted trabecular projects through edges of registered shape
% % % %%%%% This seems most likely due to the lack of accuracy in certain areas
% % % %%%%% (e.g. around ends, curvature) of the registered new tibial cortical
% % % %%%%% surface to the shape model PC scores...
% % % 
% % % %%%%% Getting somewhat better using the rigid CPD algorithm --- but perhaps
% % % %%%%% still needs work, or a better matching original tibia shape...this
% % % %%%%% also scales...but perhaps we can scale up and then scale down...
% % % 
% % % %%%%% The trick is the matching back to the original shape from the
% % % %%%%% shape model reconstructing the tibia --- the predicted trabecular
% % % %%%%% sits well inside of the 'optimised' reconstruction, but projects
% % % %%%%% through the original shape because it doesn't match the optimisation
% % % %%%%% that well...
% % % 
% % % %%%%% If including a scaling element to the rigid transformation --- then
% % % %%%%% we need to unscale back when we transform back to the original and
% % % %%%%% moving the trabecular to match this...
% % % 
% % % %%%%% Because there is such little room for error, the trabecular can still
% % % %%%%% project out when there is an OK reconstruction --- the predictive
% % % %%%%% model isn't foolproof to this...
% % %     %%%%% There's two options for me:
% % %         %%%%% (1) if the trabecular prediction always sits inside the
% % %         %%%%% optimised cortical, but not the orignal --- then there needs
% % %         %%%%% to be a transform/scale determined that alters the optimised
% % %         %%%%% surface back to the original. This transformation/scale is
% % %         %%%%% then applied to the trabecular so it sits nicely...
% % %         %%%%% (2) some sort of transformation/scale needs to be applied to
% % %         %%%%% the predicted trabecular so that it sits nicely within the
% % %         %%%%% original cortical, irrespective of what the reconstruction
% % %         %%%%% looks like...
% % % 
% % % %%% Test applying procrustes transform to new trabecular with respect to
% % % %%% registered tibia --- or rigid transform via CPD
% % % [~, proc_trabV] = procrustes(regTibV, newTrabV, 'scaling', false);



%% ...








