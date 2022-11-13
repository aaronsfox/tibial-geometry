%% This script provides a demonstration of how to generate a volumetric mesh
%  of the tibia trabecular by registering a different surface to the shape
%  model of the tibia and trabecular. The new tibial surface is matched
%  against one created by the shape model, and then the trabecular that
%  stems from these shape model component scores is created.
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

%% Options

%Flag to run all reconstructions and calculations or load in results
%Default to false as this takes time
runReconstructions = false;

%% Set-up

%Turn off warnings
warning off

%Set home directory
homeDir = pwd;

%Add supplementary code paths
addpath(genpath('supplementary'));

%Load the tibia and trabecular shape models required
load('..\..\ShapeModels\tibia-plus-trabecular\tibTrabShapeModel.mat');

%Set array to rotate views of surfaces later in reconstruction process
surfaceRot = [-90, 0, 90, 180];

%Set list to label reconstruction subplots with
viewLabel = [{'Anterior'}, {'Lateral'}, {'Posterior'}, {'Medial'}];

%Navigate to segmentation directory to get case ID names
cd('..\..\Segmentation\');

%Grab the case names based on folders
f = dir(pwd);
%Get a logical vector that tells which is a directory.
dirFlags = [f.isdir];
%Extract only those that are directories.
subFolders = f(dirFlags);
%Get only the folder names into a cell array.
caseID = {subFolders(3:end).name}; %Start at 3 to skip . and ..

%Return to home directory
cd(homeDir);

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

%%

%Check whether to run reconstructions or load results
if runReconstructions

    %Set-up waitbar for processing trabecular surfaces
    wbar_trab = waitbar(0, 'Predicting trabecular surfaces...');

    %Loop through predicted reconstructions
    %Align them in the appropriate way
    %Calulcate the mean and peak error, and Jaccard similarity
    for predictCase = 1:length(caseID)

        %Get the original surface in XYZ format
        actualV = reshape(tibTrabShapeModel.nodes(predictCase,:), ...
            [3, length(tibTrabShapeModel.mean)/3])';
        
        %Split this into the tibia and trabecular points
        %(i.e. first half of points = tibia; second half of points = trab)
        actualV_tib = actualV(1:optionStruct_tib.nb_pts,:);
        actualV_trab = actualV(optionStruct_tib.nb_pts+1:end,:);
        
        %Set-up the initial guess and bounds for the optimisation
        
        %Grab the original surface points and remove the mean
        optimData = tibTrabShapeModel.nodes(predictCase,:) - tibTrabShapeModel.mean;
        
        %Project the new data against the loadings to find the estimated scores
        %Note that these should be pretty accuracte given it's the real data
        %Set the upper and lower bounds on the score variables. Here we set them as
        %-5 / +5 standard deviations about the PC score mean (i.e. zero)        
        sdRange = 5;
        %Loop through retained PCs
        for nPC = 1:length(tibTrabShapeModel.varExplained)
            %Calculate scores and set to an initial guess variable
% % %             x0(nPC,1) = dot(tibTrabShapeModel.loadings(:,nPC), optimData);
            x0(nPC,1) = 0;
            %Set the lower bound
            lb(nPC,1) = std(tibTrabShapeModel.score(:,nPC)) * -sdRange;
            %Set the upper bound
            ub(nPC,1) = std(tibTrabShapeModel.score(:,nPC)) * sdRange;
        end
        
        %Create the function handle for fmincon
        optFunc = @(pcScores)calcReconstructionError(pcScores, tibTrabShapeModel, actualV_tib);

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
        reconstructedV_tib = optReconstructedV(1:optionStruct_tib.nb_pts,:);
        reconstructedV_trab = optReconstructedV(optionStruct_tib.nb_pts+1:end,:);

% % %         %Visualise original vs. reconstructed tibia
% % %         cFigure; hold on;
% % %         gpatch(tibTrabShapeModel.F1, actualV_tib, 'gw', 'none', 0.3);
% % %         gpatch(tibTrabShapeModel.F1, reconstructedV_tib, 'rw', 'k', 1);
% % % % % %         plotV(reconstructedV_trab, 'b.');
% % %         axisGeom; camlight headlight
        
        %Examine the number of trabecular points inside the tibia with the
        %optimised reconstructed surfaces
        trabPtsIn = intriangulation(actualV_tib, tibTrabShapeModel.F1, reconstructedV_trab, 0);

        %Extract the initial number of points outside        
        trabPointsOutsideInitial(predictCase,1) = sum(~trabPtsIn);
        
        %Rescale trabecular if not all points are inside tibia
        if sum(~trabPtsIn) > 0
            reconstructedV_trab = alignTrabecularModel(actualV_tib, reconstructedV_trab, tibTrabShapeModel, trabPtsIn);
        end
        
        %Determine trabecular points inside vs. outside tibia after alignment
        trabPtsInNew = intriangulation(actualV_tib, tibTrabShapeModel.F1, reconstructedV_trab, 0);

        %Extract the initial number of points outside        
        trabPointsOutsideAfter(predictCase,1) = sum(~trabPtsInNew);

        %Calculate the Jaccard index between actual and predicted trabecular
        jaccardSimilarity(predictCase,1) = calcJaccardTrabecular(actualV_trab, reconstructedV_trab, tibTrabShapeModel);

        %Calculate point error distance
        pointErrorDist(predictCase,:) = distancePoints3d(actualV_trab, reconstructedV_trab);

        %Calculate the mean error
        pointErrorDistMean(predictCase,1) = mean(pointErrorDist(predictCase,:));

        %Calculate the peak error
        pointErrorDistMax(predictCase,1) = max(pointErrorDist(predictCase,:));

        %Convert distance error to colour scales for visualisation
        pointErrorDistColF(predictCase,:) = vertexToFaceMeasure(tibTrabShapeModel.F2, pointErrorDist(predictCase,:)');
        pointErrorDistColV(predictCase,:) = faceToVertexMeasure(tibTrabShapeModel.F2, actualV_trab, pointErrorDistColF(predictCase,:)');
        
        %Create visualisation of error
        %Use subplots to create different perspectives
        cFigure; hold on;
        %Loop through four views to create subplot
        for viewNo = 1:4    
            %Create subplot for current view
            subplot(1,4,viewNo);
            %Add greyed out original surface
            hpOrig = gpatch(tibTrabShapeModel.F2, actualV_trab, [200/255 200/255 200/255], 'none', 0.5);
            %Add colormapped reconstruction
            hpPred = gpatch(tibTrabShapeModel.F2, reconstructedV_trab, pointErrorDistColV(predictCase,:)', 'none', 1);
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
        export_fig(['figures\predictedReconstructionError\',caseID{predictCase},'_predictedReconstructionErrorMap.png'],'-m2');
        close

        %Save predicted surface as an STL
        %%%% IMPORTANT NOTE: THESE STL FILES MUST BE CONSIDERED 'NEW' SURFACES
        %%%% AS THE POINTS GET REORDERED AND NO LONGER ALIGN WITH THE ORIGINALS
        stlStruct.solidNames = {['predicted-trabecular-',caseID{predictCase}]}; %names of parts
        stlStruct.solidVertices = {reconstructedV_trab}; %Vertices
        stlStruct.solidFaces = {tibTrabShapeModel.F2}; %Faces
        stlStruct.solidNormals={[]};

        %Export STLs
        export_STL_txt(['predictedSurfaces\',caseID{predictCase},'-predicted-trabecular.stl'], stlStruct);

        %Store data in array
        predictedSurfaces.V(:,:,predictCase) = reconstructedV_trab;
        predictedSurfaces.F(:,:,predictCase) = tibTrabShapeModel.F2;

        %Update waitbar
        waitbar(predictCase/length(caseID), wbar_trab);

    end

    %Close waitbar
    close(wbar_trab);

    %Place calculations into structure and tables
    trabecularErrorSummary.errorSummaryTable = table(pointErrorDistMean, ...
        pointErrorDistMax, ...
        jaccardSimilarity, ...
        'VariableNames', {'pointErrorDistMean', 'pointErrorDistMax', 'jaccardSimilarity'}, ...
        'RowNames', caseID);

    %Export table to file
    writetable(trabecularErrorSummary.errorSummaryTable, 'results\errorSummaryTable.csv', ...
        'WriteRowNames', true);
    
    %Place results in structure
    trabecularErrorSummary.trabPointsOutsideInitial = trabPointsOutsideInitial;
    trabecularErrorSummary.trabPointsOutsideAfter = trabPointsOutsideAfter;
    trabecularErrorSummary.jaccardSimilarity = jaccardSimilarity;
    trabecularErrorSummary.pointErrorDist = pointErrorDist;
    trabecularErrorSummary.pointErrorDistMean = pointErrorDistMean;
    trabecularErrorSummary.pointErrorDistMax = pointErrorDistMax;
    trabecularErrorSummary.pointErrorDistColF = pointErrorDistColF;
    trabecularErrorSummary.pointErrorDistColV = pointErrorDistColV;
    trabecularErrorSummary.predictedSurfaces.V = predictedSurfaces.V;
    trabecularErrorSummary.predictedSurfaces.F = predictedSurfaces.F;
    
    %Save error summary as mat
    save('results\trabecularErrorSummary.mat', 'trabecularErrorSummary');

else
    
    %Load pre-processed results
    load('results\trabecularErrorSummary.mat');
    
end

%Calculate and display mean and 95% CI's for the mean and peak error
%Mean & SD for mean error
meanError_m = mean(trabecularErrorSummary.pointErrorDistMean);
meanError_sd = std(trabecularErrorSummary.pointErrorDistMean);
meanError_lower95 = meanError_m - (1.96 * (meanError_sd / sqrt(length(caseID))));
meanError_upper95 = meanError_m + (1.96 * (meanError_sd / sqrt(length(caseID))));
%Mean & SD for maxerror
maxError_m = mean(trabecularErrorSummary.pointErrorDistMax);
maxError_sd = std(trabecularErrorSummary.pointErrorDistMax);
maxError_lower95 = maxError_m - (1.96 * (maxError_sd / sqrt(length(caseID))));
maxError_upper95 = maxError_m + (1.96 * (maxError_sd / sqrt(length(caseID))));
%Jaccard similarity
jaccard_m = mean(trabecularErrorSummary.jaccardSimilarity);
jaccard_sd = std(trabecularErrorSummary.jaccardSimilarity);
jaccard_lower95 = jaccard_m - (1.96 * (jaccard_sd / sqrt(length(caseID))));
jaccard_upper95 = jaccard_m + (1.96 * (jaccard_sd / sqrt(length(caseID))));
%Display
disp(['Reconstruction mean point error (mean +/- 95% CIs) = ',num2str(round(meanError_m,2)), ...
    '[',num2str(round(meanError_lower95,2)),',',num2str(round(meanError_upper95,2)),']']);
disp(['Reconstruction max point error (mean +/- 95% CIs) = ',num2str(round(maxError_m,2)), ...
    '[',num2str(round(maxError_lower95,2)),',',num2str(round(maxError_upper95,2)),']']);
disp(['Reconstruction Jaccard Index (mean +/- 95% CIs) = ',num2str(round(jaccard_m,3)), ...
    '[',num2str(round(jaccard_lower95,3)),',',num2str(round(jaccard_upper95,3)),']']);

%Create summary figure of error data
%Create figure
hfErr = cFigure; hold on
hfErr.Units = 'centimeters';
hfErr.Position = [5, 5, 24, 6];
%Plot mean error data
%Create subplot
subplot(1,3,1); hold on
%Plot mean as dashed line & CIs as dotted lines
plot([1,length(caseID)], [meanError_m, meanError_m], ...
    'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
plot([1,length(caseID)], [meanError_lower95, meanError_lower95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
plot([1,length(caseID)], [meanError_upper95, meanError_upper95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
%Plot participant data
scatter(linspace(1,length(caseID), length(caseID)), trabecularErrorSummary.pointErrorDistMean, ...
    15, 'k', 'filled');
%Set labels
xlabel('Participants'); ylabel('Error (mm)'); title('Mean Position Error');
%Remove ticks
ax = gca(); ax.XTick = [];
%Plot max error data
%Create subplot
subplot(1,3,2); hold on
%Plot mean as dashed line & CIs as dotted lines
plot([1,length(caseID)], [maxError_m, maxError_m], ...
    'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
plot([1,length(caseID)], [maxError_lower95, maxError_lower95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
plot([1,length(caseID)], [maxError_upper95, maxError_upper95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
%Plot participant data
scatter(linspace(1,length(caseID), length(caseID)), trabecularErrorSummary.pointErrorDistMax, ...
    15, 'k', 'filled');
%Set labels
xlabel('Participants'); ylabel('Error (mm)'); title('Peak Position Error');
%Remove ticks
ax = gca(); ax.XTick = [];
%Plot Jaccard similarity data
%Create subplot
subplot(1,3,3); hold on
%Plot mean as dashed line & CIs as dotted lines
plot([1,length(caseID)], [jaccard_m, jaccard_m], ...
    'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
plot([1,length(caseID)], [jaccard_lower95, jaccard_lower95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
plot([1,length(caseID)], [jaccard_upper95, jaccard_upper95], ...
    'Color', 'k', 'LineWidth', 1, 'LineStyle', ':');
%Plot participant data
scatter(linspace(1,length(caseID), length(caseID)), trabecularErrorSummary.jaccardSimilarity, ...
    15, 'k', 'filled');
%Set labels
xlabel('Participants'); ylabel('Jaccard Index (0-1)'); title('Jaccard Similarity');
%Remove ticks
ax = gca(); ax.XTick = [];
%Export summary figure
export_fig('figures\errorSummary\errorSummaryFigure.png','-m5');
close(hfErr);

%% Apply trabecular prediction to a new surface
%  Note that here we randomly select one surface from the dataset, but this
%  can be changed to another from the dataset. Alternatively, any .stl
%  surface similar to this example could be used.

%Select a random surface from the Nolte et al. dataset
rng(12345)
sampleNo = randi([1,35],1,1);

%Set sample ID based on random selection
if sampleNo < 10
    pID = ['C0',num2str(sampleNo),'RTF'];
else
    pID = ['C',num2str(sampleNo),'RTF'];
end

%Load data
%Tibia-fibula
[tibFibSTLstruct] = import_STL(['sampleSurfaces\Nolte2016\',pID,'.stl']);
tibFibF = tibFibSTLstruct.solidFaces{1}; %Faces
tibFibV = tibFibSTLstruct.solidVertices{1}; %Vertices
[tibFibF,tibFibV] = mergeVertices(tibFibF,tibFibV);

%Split the surfaces to grab the tibia only

%Group the vertices and faces from the combined shape model
[groupIndexVertices,groupIndexFaces] = groupVertices(tibFibF, tibFibV,0);

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

%Remesh the surfaces to match the shape model points
[tibiaF,tibiaV] = ggremesh(tibiaF, tibiaV, optionStruct_tib);

%Generate the predicted trabecular points from shape models
[newPredictedTrabV, newTibiaF, newTibiaV] = ...
    createTrabecularFromSample(tibiaF, tibiaV, tibTrabShapeModel, optionStruct_tib.nb_pts);

%Review newly fixed trabecular within original tibia
cFigure; hold on;
%Loop through four views to create subplot
for viewNo = 1:4    
    %Create subplot for current view
    subplot(1,4,viewNo);
    %Add greyed out tibia surface
    hpTib = gpatch(newTibiaF, newTibiaV, 'kw', 'none', 0.3);
    %Add trabecular surface
    hpTrab = gpatch(tibTrabShapeModel.F2, newPredictedTrabV, 'rw', 'k', 1);
    %Set axis view
    axis equal; axis tight; view(0,90);
    rotate(hpTib,[0 1 0], surfaceRot(viewNo));
    rotate(hpTrab,[0 1 0], surfaceRot(viewNo));
    %Set axis parameters
    camlight headlight; axis off
    %Add title
    title([viewLabel{viewNo},' View'], 'FontSize', 12);
end
% % % %Add figure title
% % % sgtitle('Predicted trabecular surface (red) from new tibia sample (grey)',...
% % %         'FontSize', 16, 'FontWeight', 'bold')

%Export figure
export_fig(['figures\newTrabecularPredictions\',pID,'_predictedTrabecular.png'],'-m2');
close

%%%%% NOTE: a final step here might be to transform the predicted
%%%%% trabecular back to the original surfaces coordinates via a rigid
%%%%% transformation so the trabecular fits the position of the original
%%%%% surface.

%% ----- end of generateTrabecularFromSurface.m ----- %%