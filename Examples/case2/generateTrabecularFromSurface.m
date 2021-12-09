%% This script provides a demonstration of how to generate a volumetric mesh
%  of the tibia trabecular by registering a different surface to the shape
%  model and then combining this with our regresiion model for the
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

%Set reference case number
refCase = '147211';

%Reorder case list to make target mesh first
caseID = [caseID(find(contains(caseID, refCase))), ...
    caseID(1:find(contains(caseID, refCase))-1), ...
    caseID(find(contains(caseID, refCase))+1:end)];

%Return to home directory
cd(homeDir);

%Set options for remeshing
optionStruct_tib.nb_pts = 3500; %Set desired number of points
optionStruct_tib.disp_on = 0; % Turn off command window text display
% % % optionStruct_fib.nb_pts = 2000; %Set desired number of points
% % % optionStruct_fib.disp_on = 0; % Turn off command window text display

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

%Loop through predicted reconstructions, calculate the mean and peak error,
%and visualise the reconstruction accuracy
for predictCase = 1:length(leaveOut)
    
    %Get the original surface in XYZ format
    actualV = reshape(trabShapeModel.nodes(predictCase,:), ...
        [3, length(trabShapeModel.mean)/3])';
    
    %Get the predicted reconstruction in XYZ format
    predictedV = reshape(predictReconstructed(predictCase,:) + trabShapeModel.mean, ...
        [3, length(trabShapeModel.mean)/3])';
    
    %Calculate point error distance
    pointErrorDist(predictCase,:) = distancePoints3d(actualV, predictedV);

    %Calculate the mean error
    pointErrorDistMean(predictCase,1) = mean(pointErrorDist(predictCase,:));
    
    %Calculate the peak error
    pointErrorDistMax(predictCase,1) = max(pointErrorDist(predictCase,:));
    
    %Convert distance error to colour scales for visualisation
    pointErrorDistColF(predictCase,:) = vertexToFaceMeasure(trabShapeModel.F, pointErrorDist(predictCase,:)');
    pointErrorDistColV(predictCase,:) = faceToVertexMeasure(trabShapeModel.F, actualV, pointErrorDistColF(predictCase,:)');
    
    %Create visualisation of error
    %Use subplots to create different perspectives
    cFigure; hold on;
    %Loop through four views to create subplot
    for viewNo = 1:4    
        %Create subplot for current view
        subplot(1,4,viewNo);
        %Add greyed out original surface
        hpOrig = gpatch(trabShapeModel.F, actualV, [200/255 200/255 200/255], 'none', 0.5);
        %Add colormapped reconstruction
        hpPred = gpatch(trabShapeModel.F, predictedV, pointErrorDistColV(predictCase,:)', 'none', 1);
        %Interpolate colouring for smoothness
        hpPred.FaceColor = 'Interp'; colormap viridis
        %Set axis view
        axis equal; axis tight; view(0,90);
        rotate(hpOrig,[0 1 0], surfaceRot(viewNo));
        rotate(hpPred,[0 1 0], surfaceRot(viewNo));
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
    export_fig(['figures\predictedReconstructionError\',caseID{predictCase},'_predictedReconstructionErrorMap.png'],'-m1');
	close
    
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

%% Apply trabecular prediction to new surface

%Load the tibia and fibula surface meshes
%Tibia
[tibiaSTLstruct] = import_STL('sampleSurfaces\rightTibia.stl');
tibiaF = tibiaSTLstruct.solidFaces{1}; %Faces
tibiaV = tibiaSTLstruct.solidVertices{1}; %Vertices
[tibiaF,tibiaV] = mergeVertices(tibiaF,tibiaV);
% % % %Fibula
% % % [fibulaSTLstruct] = import_STL('sampleSurfaces\rightFibula.stl');
% % % fibulaF = fibulaSTLstruct.solidFaces{1}; %Faces
% % % fibulaV = fibulaSTLstruct.solidVertices{1}; %Vertices
% % % [fibulaF,fibulaV] = mergeVertices(fibulaF,fibulaV);

%Remesh the surfaces to match the shape model points
[tibiaF,tibiaV] = ggremesh(tibiaF, tibiaV, optionStruct_tib);

% % % %Visualise
% % % cFigure; hold on
% % % gpatch(tibiaF,tibiaV,'gw','k');
% % % % % % gpatch(fibulaF,fibulaV,'bw','k');
% % % axisGeom; camlight headlight
% % % title('Imported Tibia Surface');

%%%%% Don't think we need fibula???


%% Update below...

%%%% TODO: register using the CPD algorithm...
%%%% TODO: create optimisation to minimise point error by manipulating
%%%% shape model scores...

%% Align new surface to tibia shape model

%Remesh the new surface to match the points of the shape model mean
%Set options
optionStruct_tib.nb_pts = length(tibiaShapeModel.meanPts); %Set desired number of points
optionStruct_tib.disp_on = 0; % Turn off command window text display
%Remesh
[tibiaF,tibiaV] = ggremesh(tibiaF, tibiaV, optionStruct_tib);

%It helps the algorithm to get the surfaces roughly aligned, so perform
%some initial transformations to get these objects in a similar orientation
%and position. The rotations and transformations applied here required a
%bit of estimation/testing.

%Rotate -90 degrees about the x-axis
%Create rotation transform
rotX = createRotationOx(deg2rad(-90));
%Apply to points
for ptNo = 1:length(tibiaV)
    tibiaV(ptNo,:) = transformPoint3d(tibiaV(ptNo,:),rotX);    
end

%Rotate 90 degrees about the y-axis
%Create rotation transform
rotY = createRotationOy(deg2rad(90));
%Apply to points
for ptNo = 1:length(tibiaV)
    tibiaV(ptNo,:) = transformPoint3d(tibiaV(ptNo,:),rotY);    
end

%Identify the magnitude required for a rough translation by comparing the
%means of the points
%Calculate mean difference
surfacePtMean = mean(tibiaV);
meanPtMean = mean(tibiaShapeModel.meanPts);
transDist = meanPtMean - surfacePtMean;
%Create translation matrix
transMat = createTranslation3d(transDist);

%Apply translation to points
for ptNo = 1:length(tibiaV)
    tibiaV(ptNo,:) = transformPoint3d(tibiaV(ptNo,:),transMat);    
end

%To map the new surface to the shape model it's important to identify the
%corresponding points on the surface to that of the shape model.

%Perform non-rigid registration using the ICP algorithm to assist in
%identifying corresponding points. This also includes an initial step to
%rigidly align the surfaces.
[registered,targetV,targetF] = nonrigidICPv1(meanV, tibiaV, ...
    meanF, tibiaF, 10, 0);

% % % %Visualise mean, original and registered meshes
% % % cFigure; hold on;
% % % gpatch(meanF, meanV, 'bw', 'k', 0.3);
% % % gpatch(tibiaF, tibiaV, 'rw', 'k', 0.3);
% % % gpatch(tibiaF, registered, 'gw', 'k', 0.3);
% % % axisGeom;

%Identify the matching points in the shape model mean against the
%registration to identify the corresponding point indices
regSortIdx = knnsearch(registered, tibiaShapeModel.meanPts);

%Sort the registered points so they align with the shape model mean points
registeredSorted = registered(regSortIdx,:);

%Sort the original surface points so they align with the shape model
tibiaV_shapeModel = tibiaV(regSortIdx,:);

% % % %Visualise to confirm all points align with mean
% % % cFigure; hold on;
% % % gpatch(meanF, tibiaShapeModel.meanPts(ptSortIdx,:), 'bw', 'k', 0.3);
% % % gpatch(meanF, tibiaV_shapeModel(ptSortIdx,:), 'rw', 'k', 0.3);
% % % gpatch(meanF, registeredSorted(ptSortIdx,:), 'gw', 'k', 0.3);
% % % axisGeom;

%At this point the points in the tibiaV_shapeModel correspond to the
%reshaped mean points in the shape model dataset. Given this, we can now
%map the new surface points to the shape model to extract the model weights
%that would generate this surface

%% Extract shape model weights for the tibia surface

%Project the new data points on to the modes of the shape model

%Reshape the surface points to a single column array
newData = reshape(tibiaV_shapeModel', [1, length(tibiaV_shapeModel)*3]);

%Subtract the mean from the new data
newData_mRemoved = newData - tibiaShapeModel.mean;

%Set the number of PCs to use. Here we define a cut-off for the tibia shape
%model at a cumulative variance > 97% to match what is used in the
%regression model
for nPC = 1:length(tibiaShapeModel.weights)
    tibiaShapeModel.varExplained(nPC) = tibiaShapeModel.weights(nPC) / sum(tibiaShapeModel.weights);
end
usePCs = find(cumsum(tibiaShapeModel.varExplained) > 0.97, 1);

%Get the feature vector from the shape model
f = tibiaShapeModel.modes(:, 1:usePCs);

%Project to find the weights
for nPC = 1:usePCs
    %Calculate weight
    w(nPC,1) = dot(f(:,nPC), newData_mRemoved);
    %Convert to normalised weight format
    wNorm(nPC,1) = (w(nPC,1) - mean(tibiaShapeModel.projectedWeights(nPC,:))) / ...
        std(tibiaShapeModel.projectedWeights(nPC,:));       
end

%Confirm that the extracted weights for new surface are appropriate by
%reconstructing and comparing it to the original surface

%Reconstruct the tibia points for new surface
for nPC = 1:length(wNorm)
    
    %Get the current PC weight
    reconW = wNorm(nPC) * sqrt(tibiaShapeModel.weights(nPC));
    
    %Get the reconstruction points for current PC
    shapePts(nPC,:) = reconW * tibiaShapeModel.modes(:,nPC)';
    
end

%Sum the calculated component points and add to the mean shape
reconSurface = sum(shapePts,1) + tibiaShapeModel.mean;

%Reshape points to visualise
reconSurfaceV = transpose(reshape(reconSurface, [3, length(reconSurface)/3]));

%Visualise original vs. reconstructed
cFigure; hold on;
gpatch(meanF, tibiaV_shapeModel(ptSortIdx,:), 'bw', 'k', 0.3);
gpatch(meanF, reconSurfaceV(ptSortIdx,:), 'rw', 'k', 0.3);
axisGeom;

%%%%% TODO: points are a little off - does the orientation and positioning
%%%%% of the new surface need to first be aligned in the same way based on
%%%%% the landmarks we used for the shape model? --- better alignment in
%%%%% the first place like with a rigid registration or something??? There
%%%%% was no real alignment in the first place...

%%%%% Nonetheless, the reconstruction seems OK - the positioning of the
%%%%% surface is not ideal though...even with rigid alignment, it still
%%%%% struggles a little bit - but perhaps this is simply the error in
%%%%% application to new data

%Calculate error between reconstructed and original surface
reconError = distancePoints3d(tibiaV_shapeModel(ptSortIdx,:), reconSurfaceV(ptSortIdx,:));

%Map errors on surface
%Convert distances to face measures
reconErrorC = vertexToFaceMeasure(meanF, reconError);
%Visualise
cFigure; hold on;
gpatch(meanF, reconSurfaceV(ptSortIdx,:), reconErrorC, 'k');
colormap viridis; colorbar
axisGeom;

%%%%% I don't necessarily think the errors and slightly wonky shape will be
%%%%% a problem IF we can use the original tibia surface with the predicted
%%%%% trabecular --- i.e. if we can align the predicted trabecular shape
%%%%% with the original surface, then everything is fine, right???
%%%%%
%%%%% To do this, we subsequently need a good approach for re-alligning the
%%%%% predicted trabecular to the original surface --- perhaps something to
%%%%% do with a rigid transform or use the two mean shape model surfaces as
%%%%% an appropriate use. If this is the case, then we probably don't need
%%%%% to use the tibia-fibula model at all, considering it probably won't
%%%%% actually help prediction of a trabecular surface that well...


%%





%Set options for remeshing

optionStruct_fib.nb_pts = 2500; %Set desired number of points
optionStruct_fib.disp_on = 0; % Turn off command window text display













