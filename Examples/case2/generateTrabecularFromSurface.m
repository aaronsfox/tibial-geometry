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

% % % %Create a scatter subplot of the predicted trabecular PCs
% % % f = figure; f.Position = [100 100 1080 800];
% % % for predictPC = 1:trabShapeModel.retainPCs
% % %     
% % %     %Create subplot
% % %     subplot(2,2,predictPC); hold on
% % %     
% % %     %Plot data
% % %     scatter(trabShapeModel.score(:,predictPC), predictTrabPCs(:,predictPC), ...
% % %         'k', 'filled');
% % %     
% % %     %Plot line of equality
% % %     ax = gca;
% % %     plot([ax.XLim(1), ax.XLim(2)], [ax.XLim(1), ax.XLim(2)], 'r--');
% % %     
% % %     %Add title
% % %     title(['Trabecular Shape Model PC',num2str(predictPC)]);
% % %     
% % %     %Add labels
% % %     xlabel('Actual Shape Model Score');
% % %     ylabel('Predicted Shape Model Score');
% % %     
% % % end

%Reconstruct each case using the predicted PCs
predictReconstructed = predictTrabPCs(:,1:trabShapeModel.retainPCs) * ...
    trabShapeModel.loadings(:,1:trabShapeModel.retainPCs)';

%Check whether to run reconstructions or load results
if runReconstructions

    %Set-up waitbar for processing trabecular surfaces
    wbar_trab = waitbar(0, 'Predicting trabecular surfaces...');

    %Loop through predicted reconstructions
    %Align them in the appropriate way
    %Calulcate the mean and peak error, and Jaccard similarity
    for predictCase = 1:length(leaveOut)

        %Get the original surface in XYZ format
        actualV = reshape(trabShapeModel.nodes(predictCase,:), ...
            [3, length(trabShapeModel.mean)/3])';

        %Get the predicted reconstruction in XYZ format
        predictedV = reshape(predictReconstructed(predictCase,:) + trabShapeModel.mean, ...
            [3, length(trabShapeModel.mean)/3])';

        %Get the tibial surface for alignment
        tibiaV = reshape(tibiaShapeModel.nodes(predictCase,:), ...
            [3, length(tibiaShapeModel.mean)/3])';

        %Check the need to scale trabecular to fit inside tibia

        %Determine trabecular points inside vs. outside tibia
        trabPtsIn = intriangulation(tibiaV, tibiaShapeModel.F, predictedV, 0);

        %Extract the initial number of points outside        
        trabPointsOutsideInitial(predictCase,1) = sum(~trabPtsIn);

        %Rescale trabecular if not all points are inside tibia
        if sum(~trabPtsIn) > 0
            predictedV = alignTrabecularModel(tibiaV, tibiaShapeModel, predictedV, trabShapeModel, trabPtsIn);
        end

        %Determine trabecular points inside vs. outside tibia after alignment
        trabPtsInNew = intriangulation(tibiaV, tibiaShapeModel.F, predictedV, 0);

        %Extract the initial number of points outside        
        trabPointsOutsideAfter(predictCase,1) = sum(~trabPtsInNew);

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

        %Save predicted surface as an STL
        %%%% IMPORTANT NOTE: THESE STL FILES MUST BE CONSIDERED 'NEW' SURFACES
        %%%% AS THE POINTS GET REORDERED AND NO LONGER ALIGN WITH THE ORIGINALS
        stlStruct.solidNames = {['predicted-trabecular-',caseID{predictCase}]}; %names of parts
        stlStruct.solidVertices = {predictedV}; %Vertices
        stlStruct.solidFaces = {trabShapeModel.F}; %Faces
        stlStruct.solidNormals={[]};

        %Export STLs
        export_STL_txt(['predictedSurfaces\',caseID{predictCase},'-predicted-trabecular.stl'], stlStruct);

        %Store data in array
        predictedSurfaces.V(:,:,predictCase) = predictedV;
        predictedSurfaces.F(:,:,predictCase) = trabShapeModel.F;

        %Update waitbar
        waitbar(predictCase/length(leaveOut), wbar_trab);

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
meanError_lower95 = meanError_m - (1.96 * (meanError_sd / sqrt(length(leaveOut))));
meanError_upper95 = meanError_m + (1.96 * (meanError_sd / sqrt(length(leaveOut))));
%Mean & SD for maxerror
maxError_m = mean(trabecularErrorSummary.pointErrorDistMax);
maxError_sd = std(trabecularErrorSummary.pointErrorDistMax);
maxError_lower95 = maxError_m - (1.96 * (maxError_sd / sqrt(length(leaveOut))));
maxError_upper95 = maxError_m + (1.96 * (maxError_sd / sqrt(length(leaveOut))));
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
export_fig('figures\errorSummary\errorSummaryFigure.png','-m2');
close(hfErr);

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
    
    %Get linear model as 
    
    %Export model coefficients to table
    coeffTable = finalLinearModel{predictPC}.Coefficients;
    coeffTable.Properties.RowNames = {'Intercept';
        'Tibia Shape Model PC1';
        'Tibia Shape Model PC2';
        'Tibia Shape Model PC3';
        'Tibia Shape Model PC4';
        'Tibia Shape Model PC5'};
    coeffTable.Properties.VariableNames = {'Estimate', 'Standard Error', 't-statistic', 'p-value'};
    writetable(coeffTable, ['results\linearModel_trabecularPC',num2str(predictPC),'.csv'], ...
        'WriteRowNames', true);
    
    %Print summary diagnostics to file
    fid = fopen(['results\linearModelResults_trabecularPC',num2str(predictPC),'.txt'], 'wt');
    T = evalc('disp(finalLinearModel{predictPC})');
    fprintf(fid, T);
    fclose(fid);
    
% % %     %Display model output
% % %     finalLinearModel{predictPC}

% % %     %Plot model
% % %     figure; plot(finalLinearModel{predictPC});

end

%% Apply trabecular prediction to a new surface

%Select a random surface from the Nolte et al. dataset
%%%% NOTE: this can be changed to select a different surface
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

% % % %Visualise
% % % cFigure; hold on
% % % gpatch(tibiaF,tibiaV,'gw', 'none')
% % % gpatch(fibulaF,fibulaV,'bw', 'none')
% % % axisGeom; camlight headlight

%Remesh the surfaces to match the shape model points
[tibiaF,tibiaV] = ggremesh(tibiaF, tibiaV, optionStruct_tib);

%Generate the predicted trabecular points from shape models
[newPredictedTrabV, newTibiaF, newTibiaV] = ...
    createTrabecularFromSample(tibiaF, tibiaV, tibiaShapeModel, trabShapeModel, finalLinearModel);

%Review newly fixed trabecular within original tibia
cFigure; hold on;
%Loop through four views to create subplot
for viewNo = 1:4    
    %Create subplot for current view
    subplot(1,4,viewNo);
    %Add greyed out tibia surface
    hpTib = gpatch(newTibiaF, newTibiaV, 'kw', 'none', 0.3);
    %Add trabecular surface
    hpTrab = gpatch(trabShapeModel.F, newPredictedTrabV, 'rw', 'k', 1);
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
export_fig(['figures\newTrabecularPredictions\',pID,'_predictedTrabecular.png'],'-m1');
close

%%%%% NOTE: a final step here might be to transform the predicted
%%%%% trabecular back to the original surfaces coordinates via a rigid
%%%%% transformation

%% ----- end of generateTrabecularFromSurface.m ----- %%