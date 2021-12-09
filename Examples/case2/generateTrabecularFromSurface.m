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

%Add supplementary code paths
addpath(genpath('supplementary'));

%Load the tibia and fibula surface meshes
%Tibia
[tibiaSTLstruct] = import_STL('rightTibia.stl');
tibiaF = tibiaSTLstruct.solidFaces{1}; %Faces
tibiaV = tibiaSTLstruct.solidVertices{1}; %Vertices
[tibiaF,tibiaV] = mergeVertices(tibiaF,tibiaV);
%Fibula
[fibulaSTLstruct] = import_STL('rightFibula.stl');
fibulaF = fibulaSTLstruct.solidFaces{1}; %Faces
fibulaV = fibulaSTLstruct.solidVertices{1}; %Vertices
[fibulaF,fibulaV] = mergeVertices(fibulaF,fibulaV);

%% Create regression model between tibia and trabecular

%Regression model for trabecular PC

%%%% TODO: cross validation k-fold approach (70:30 split) to train
%%%% regression model???

%Create prediction matrix
%Use retained PCs as test
X = tibiaShapeModel.score(:,1:tibiaShapeModel.retainPCs);
% % % X = zscore(tibiaShapeModel.score(:,1:tibiaShapeModel.retainPCs));

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


%% Mesh processing

%Remesh the surfaces to match the shape model points

% % % [fibulaF,fibulaV] = ggremesh(fibulaF, fibulaV, optionStruct_fib);

% % % %Visualise
% % % cFigure; hold on
% % % gpatch(tibiaF,tibiaV,'gw','k');
% % % gpatch(fibulaF,fibulaV,'bw','k');
% % % axisGeom; camlight headlight
% % % title('Imported Tibia-Fibula Surfaces');

%% Tibia only

%% Generate regression model between tibia and trabecular shape models

%Load tibia shape model
tibiaShapeModel = load('..\..\ShapeModels\tibia\shape_model\tibia.pc.mat');

%Reshape mean
tibiaShapeModel.meanPts = transpose(reshape(tibiaShapeModel.mean,...
    [3, length(tibiaShapeModel.mean)/3]));

%Load mean STL
[stlStruct] = import_STL('..\..\ShapeModels\tibia\shape_model\tibia_mean.stl');
meanF = stlStruct.solidFaces{1}; %Faces
meanV = stlStruct.solidVertices{1}; %Vertices
[meanF,meanV] = mergeVertices(meanF,meanV);

%Create point sorting index to match up to faces
ptSortIdx = knnsearch(tibiaShapeModel.meanPts, meanV);

%Normalise each row of projected weights to get them in normalised/SD form
for nPC = 1:length(tibiaShapeModel.weights)
    tibiaShapeModel.normProjectedWeights(nPC,:) = zscore(tibiaShapeModel.projectedWeights(nPC,:));
end

%%%%% TODO: add regression model approachh here...

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













