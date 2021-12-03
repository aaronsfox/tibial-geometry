

%%%%% NOTE: fibula shape model looking better now --- trabecular still
%%%%% potentially has an issue, may now be due to the registration mapping
%%%%% back to the first case though. The first case has some iffy areas
%%%%% around the bottom of the surface, and this results in the fitted
%%%%% meshes looking a bit wonky at that area, and the shape model smaller
%%%%% components overly focusing on this. Hence, if we register the
%%%%% surfaces to a more 'normal' looking trabecular surface --- it may fix
%%%%% this issue...?
%%%%%
%%%%% Trabecular shape model variations appear much better when using an
%%%%% appropriate surface to register the others too...nicely done!
%%%%%
%%%%% Registering tibia alone to a different one dramatically changed %
%%%%% variance for initial components??? Same for tibia-fibula...
    %%%%% The trabecular seems better, but it may just be better than the
    %%%%% 102480 example --- as the tibia and tibia-fibula models aren't
    %%%%% that great...
    
%%%%% Tibia and tibia-fibula seem OK when using the alphanumeric ordered
%%%%% set --- trabecular still causing issues though...?


cd('..\ShapeModels\tibia\shape_model');
% cd('..\ShapeModels\trabecular\shape_model');

%Load matlab version
load('tibia.pc.mat');
% load('trabecular.pc.mat');

%Rename mean variable to avoid function clash
meanShape = mean;
clear mean

%Reshape mean
meanPts = transpose(reshape(meanShape, [3, length(meanShape)/3]));

% % % %Sort by X-axis
% % % %Stash away the indices of this sorting process so that any newly created
% % % %points can be reshaped and ordered in the same way to match up with the
% % % %faces from the mean STL file created (i.e. the STL file sorts the vertices
% % % %based on the X-axis values in ascending order)
% % % [meanPts, ptSortIdx] = sortrows(meanPts);

%Replicate getWeightsBySD function
%SD = 3.0; PC = 1 (i.e. first)
plusSD = 3.0*sqrt(weights(1));
minusSD = -3.0*sqrt(weights(1));

%Replicate reconstruct function for SD
%For PC 1 @ 3SD // i.e. using the weight above
f = modes(:,1);
plus = (plusSD * f') + meanShape;
minus = (minusSD * f') + meanShape;
plusV = transpose(reshape(plus, [3, length(plus)/3]));
minusV = transpose(reshape(minus, [3, length(minus)/3]));

%Visualise points
cFigure; hold on
plotV(meanPts, 'b.');
plotV(plusV, 'g.');
plotV(minusV, 'r.');
axisGeom;

%Load mean STL
[tibiaSTLstruct] = import_STL('tibia_mean.stl');
stl_meanF = tibiaSTLstruct.solidFaces{1}; %Faces
stl_meanV = tibiaSTLstruct.solidVertices{1}; %Vertices
[stl_meanF,stl_meanV] = mergeVertices(stl_meanF,stl_meanV);

%Sorting doesn't seem 100% foolproof - so find the matching closest point
%on the mean points and store this in an indexing variable
%This ptSortIdx can therefore be used across reconstructed surface points
%after the reshape - and theoretically should align to the face indices
%coming from the imported STL
for ptNo = 1:length(meanPts)

    %Calculate current point distances to STL nodes
    distCurrPoint = distancePoints3d(stl_meanV(ptNo,:), meanPts);

    %Get the index of the min distance points
    [~, ptSortIdx(ptNo,1)] = min(distCurrPoint);

end

% % % %Sort the mean points using the indices
% % % sortMeanPts = meanPts(ptSortIdx,:);

% % % %Visualise the newly created surfaces using the existing faces
% % % cFigure; hold on
% % % %Mean
% % % subplot(1,3,1); title('Mean Surface');
% % % gpatch(stl_meanF, meanPts(ptSortIdx, :), 'bw', 'k');
% % % axisGeom;
% % % %Plus 3SD
% % % subplot(1,3,2); title('Plus 3SD');
% % % gpatch(stl_meanF, plusV(ptSortIdx, :), 'gw', 'k');
% % % axisGeom;
% % % %Minus 3SD
% % % subplot(1,3,3); title('Minus 3SD');
% % % gpatch(stl_meanF, minusV(ptSortIdx, :), 'rw', 'k');
% % % axisGeom;

%Visualise the newly created surfaces using the existing faces
cFigure; hold on
gpatch(stl_meanF, meanPts(ptSortIdx, :), 'bw', 'k', 0.3);
gpatch(stl_meanF, plusV(ptSortIdx, :), 'gw', 'k', 0.3);
gpatch(stl_meanF, minusV(ptSortIdx, :), 'rw', 'k', 0.3);
axisGeom;


%%%% Above processes seem relatively foolproof for reconstructing PCs in
%%%% Matlab from the model and then visualising using stored mean surface

%%%% Theoretically the same sort could be applied to recreating any surface
%%%% points from the model (i.e. using a participants data)...?

%% Calculate proportional variance for each PC

%Load matlab version
load('tibia.pc.mat');

%Rename mean variable to avoid function clash
meanShape = mean;
clear mean

%Print the % variance for each PC
%Get the components
for nPC = 1:length(weights)
    varExplained(nPC) = weights(nPC) / sum(weights);
    disp(['PC',num2str(nPC),': ',num2str(round(varExplained(nPC)*100,2)),'%'])
end

%Calculate and display cumulative variance
varExplainedCumulative = cumsum(varExplained);
for nPC = 1:length(weights)
    disp(['PC',num2str(nPC),': ',num2str(round(varExplainedCumulative(nPC)*100,2)),'%'])
end


%% Attempt to reconstruct case surface from principal components
%  Aim here is to test out process of creating 'new' models

%Normalise each row of projected weights to get them in normalised/SD form
for nPC = 1:length(weights)
    normProjectedWeights(nPC,:) = zscore(projectedWeights(nPC,:));
end

%Reconstruct the tibia points about the mean for participant 1
caseInd = 1;
for nPC = 1:length(weights)
    
    %Get the current PC weight
    w = normProjectedWeights(nPC, caseInd) * sqrt(weights(nPC));
    
    %Get the reconstruction points for current PC
    shapePts(nPC,:) = w * modes(:,nPC)';
    
end

%Sum the calculated component points and add to the mean shape
reconCase = sum(shapePts,1) + meanShape;

%Reshape points to visualise
reconCaseV = transpose(reshape(reconCase, [3, length(reconCase)/3]));
meanPts = transpose(reshape(meanShape, [3, length(meanShape)/3]));

% % % %Visualise points
% % % cFigure; hold on
% % % plotV(meanPts, 'b.');
% % % plotV(reconCaseV, 'r.');
% % % axisGeom;
% % % 
% % % %Visualise the newly created surfaces using the existing faces
% % % cFigure; hold on
% % % gpatch(stl_meanF, meanPts(ptSortIdx, :), 'bw', 'k', 0.3);
% % % gpatch(stl_meanF, reconCaseV(ptSortIdx, :), 'rw', 'k', 0.3);
% % % axisGeom;

%Review reconstructed points of case 1 against the mesh fitted for case

%Load the cases mesh
[caseSTLstruct] = import_STL('..\aligned_meshes\102480-tibia_rbfreg_rigidreg.stl');
stl_caseF = caseSTLstruct.solidFaces{1}; %Faces
stl_caseV = caseSTLstruct.solidVertices{1}; %Vertices
[stl_caseF,stl_caseV] = mergeVertices(stl_caseF,stl_caseV);

%Visualise comparison
% % % cFigure; hold on
% % % gpatch(stl_caseF, stl_caseV, 'bw', 'k', 0.3);
% % % gpatch(stl_meanF, reconCaseV(ptSortIdx, :), 'rw', 'k', 0.3);
% % % axisGeom;
cFigure; hold on
plotV(stl_caseV, 'b.');
plotV(reconCaseV(ptSortIdx, :), 'r.');
axisGeom;

%%%% Exact fit here - so effectively we can test to see how reducing the
%%%% number of PCs messes with the point error of the reconstructed shape
%%%% (i.e. reduce number of PCs and see what sort of error we end up with)

%Calculate differences between points of surface to reconstruction
reconstructedOrderedPts = reconCaseV(ptSortIdx, :);
pointErrors = distancePoints3d(stl_caseV, reconstructedOrderedPts);

%%%% The above code gives certain magnitudes of error despite the points
%%%% being overlayed exactly the same --- something to do with how the
%%%% points are ordered even though the faces work??? Probably because the
%%%% stl version for the case is not the same as the mean version given the
%%%% way they are probably individually sorted --- if we look for the
%%%% nearest point or even re-registered I'm sure the error will basically
%%%% be zero

%Calculate nearest point between case STL and reconstructed
for ptNo = 1:length(reconCaseV)

    %Calculate current point distances to STL nodes
    distCurrPoint = distancePoints3d(stl_caseV(ptNo,:), reconCaseV);

    %Get the index of the min distance points
    [~, caseSortIdx(ptNo,1)] = min(distCurrPoint);

end

%Calculate new distance error
pointErrors = distancePoints3d(stl_caseV, reconCaseV(caseSortIdx,:));
disp(['Average point distance: ',num2str(mean(pointErrors))]);

%% Create heat map for PC from the mean shape & a reconstruction

%Create reconstructed points for PC1 +/- 3 SD

%Load matlab version
load('tibia.pc.mat');

%Rename mean variable to avoid function clash
meanShape = mean;
clear mean

%Reshape mean
meanPts = transpose(reshape(meanShape, [3, length(meanShape)/3]));

%Replicate getWeightsBySD function
nSD = 3; PC = 1;
plusSD = nSD*sqrt(weights(PC));
minusSD = -nSD*sqrt(weights(PC));

%Replicate reconstruct function for SD
plus = (plusSD * modes(:,PC)') + meanShape;
minus = (minusSD * modes(:,PC)') + meanShape;
plusPts = transpose(reshape(plus, [3, length(plus)/3]));
minusPts = transpose(reshape(minus, [3, length(minus)/3]));

% % % %Visualise points
% % % cFigure; hold on
% % % plotV(meanPts, 'b.');
% % % plotV(plusPts, 'g.');
% % % plotV(minusPts, 'r.');
% % % axisGeom;

%Load mean STL
% [tibiaSTLstruct] = import_STL('tibia_mean.stl');
[tibiaSTLstruct] = import_STL('trabecular_mean.stl');
stl_meanF = tibiaSTLstruct.solidFaces{1}; %Faces
stl_meanV = tibiaSTLstruct.solidVertices{1}; %Vertices
[stl_meanF,stl_meanV] = mergeVertices(stl_meanF,stl_meanV);

%Sorting doesn't seem 100% foolproof - so find the matching closest point
%on the mean points and store this in an indexing variable
%This ptSortIdx can therefore be used across reconstructed surface points
%after the reshape - and theoretically should align to the face indices
%coming from the imported STL
for ptNo = 1:length(meanPts)

    %Calculate current point distances to STL nodes
    distCurrPoint = distancePoints3d(stl_meanV(ptNo,:), meanPts);

    %Get the index of the min distance points
    [~, ptSortIdx(ptNo,1)] = min(distCurrPoint);

end

%Sort the points by the index to calculate the distance between them
meanPtsSorted = meanPts(ptSortIdx, :);
plusPtsSorted = plusPts(ptSortIdx, :);
minusPtsSorted = minusPts(ptSortIdx, :);


%Calculate distances between points for heat mapping
pcPlusDist = distancePoints3d(meanPtsSorted,plusPtsSorted);
pcMinusDist = distancePoints3d(meanPtsSorted,minusPtsSorted);

%Convert distances to face measures
pcPlusC = vertexToFaceMeasure(stl_meanF,pcPlusDist);
pcMinusC = vertexToFaceMeasure(stl_meanF,pcMinusDist);

%Visualise surfaces and change in points

%Create figure
cFigure;

%Minus PC
subplot(1,4,1);
title(['PC',num2str(PC),' -',num2str(nSD),' SD'],'fontsize',10);
gpatch(stl_meanF,minusPtsSorted,'rw','none')
axisGeom; view(2); axis off
camlight headlight

%Mean
subplot(1,4,2);
title('Mean','fontsize',10);
gpatch(stl_meanF,meanPtsSorted,'bw','none')
axisGeom; view(2); axis off
camlight headlight

%Plus PC
subplot(1,4,3);
title(['PC',num2str(PC),' +',num2str(nSD),' SD'],'fontsize',10);
gpatch(stl_meanF,plusPtsSorted,'gw','none')
axisGeom; view(2); axis off
camlight headlight

%Heatmap
subplot(1,4,4);
title(['Change Map PC',num2str(PC)]);
gpatch(stl_meanF,plusPtsSorted,pcPlusC,'none');
colormap viridis; %colorbar;
axisGeom; view(2); axis off;
camlight headlight

%Reset all axes so that they are the same
%Get the limits from each axes
for aa = 1:4
    subplot(1,4,aa);
    ax = gca;
    xLim(aa,:) = ax.XLim;
    yLim(aa,:) = ax.YLim;
    zLim(aa,:) = ax.ZLim;    
end
clear aa
%Get the min and max from each
setX(1,1) = min(xLim(:,1)); setX(1,2) = max(xLim(:,2));
setY(1,1) = min(yLim(:,1)); setY(1,2) = max(yLim(:,2));
setZ(1,1) = min(zLim(:,1)); setZ(1,2) = max(zLim(:,2));
%Set the axes values
for aa = 1:4
    subplot(1,4,aa);
    ax = gca;
    ax.XLim = setX;
    ax.YLim = setY;
    ax.ZLim = setZ;
end
clear aa

%% Create animation of shape change for PC

%Load matlab version
% load('tibia.pc.mat');
% load('trabecular.pc.mat');
load('tibia-fibula.pc.mat');

%Set PC
PC = 2;

%Rename mean variable to avoid function clash
meanShape = mean;
clear mean

%Reshape mean
meanPts = transpose(reshape(meanShape, [3, length(meanShape)/3]));

%Load mean STL
% [stlStruct] = import_STL('tibia_mean.stl');
[stlStruct] = import_STL('tibia-fibula_mean.stl');
% [stlStruct] = import_STL('trabecular_mean.stl');
stl_meanF = stlStruct.solidFaces{1}; %Faces
stl_meanV = stlStruct.solidVertices{1}; %Vertices
[stl_meanF,stl_meanV] = mergeVertices(stl_meanF,stl_meanV);

%Create point sorting index to match up to faces
%%%% TODO: replace with knnsearch...
for ptNo = 1:length(meanPts)

    %Calculate current point distances to STL nodes
    distCurrPoint = distancePoints3d(stl_meanV(ptNo,:), meanPts);

    %Get the index of the min distance points
    [~, ptSortIdx(ptNo,1)] = min(distCurrPoint);

end

%Create a variable to work through 0.1 increments from -3 to +3 SD
calcSD = linspace(0.1, 3, diff([0,3])*10);

%Calculate the points at each increment
for nSD = 1:length(calcSD)
    
    %Get weights
    plusSD = calcSD(nSD)*sqrt(weights(PC));
    minusSD = -calcSD(nSD)*sqrt(weights(PC));

    %Reconstruct points
    plus = (plusSD * modes(:,PC)') + meanShape;
    minus = (minusSD * modes(:,PC)') + meanShape;
    plusPts = transpose(reshape(plus, [3, length(plus)/3]));
    minusPts = transpose(reshape(minus, [3, length(minus)/3]));
    
    %Set into structure
    posV(:,:,nSD) = plusPts(ptSortIdx,:);
    negV(:,:,nSD) = minusPts(ptSortIdx,:);
    
    %Calculate distances between points for heat mapping
    pcPlusDist = distancePoints3d(meanPts(ptSortIdx,:), plusPts(ptSortIdx,:));
    pcMinusDist = distancePoints3d(meanPts(ptSortIdx,:), minusPts(ptSortIdx,:));

    %Convert distances to face measures
    pcPlusC = vertexToFaceMeasure(stl_meanF,pcPlusDist);
    pcMinusC = vertexToFaceMeasure(stl_meanF,pcMinusDist);
    
    %Set into structure
    posC(:,nSD) = pcPlusC;
    negC(:,nSD) = pcMinusC; 
    
end

%Create the basic view to store graphics and initiate animation
hf = cFigure;
%Positive SD
subplot(1,2,1);
title(['Plus SD for PC',num2str(PC)]);
hp_pos = gpatch(stl_meanF, meanPts(ptSortIdx,:),posC(:,1),'none',1);
axisGeom; camlight headlight
axis(axisLim([posV;negV])); %set axis limits based on values
axis off
%Negative SD
subplot(1,2,2);
title(['Minus SD for PC',num2str(PC)]);
hp_neg = gpatch(stl_meanF, meanPts(ptSortIdx,:),negC(:,1),'none',1);
axisGeom; camlight headlight
axis(axisLim([posV;negV])); %set axis limits based on values
axis off

%Set up animation
animStruct.Time = calcSD;
%Loop through SD points
for nSD = 1:length(calcSD)
    
    %Set entries into animation structure
    animStruct.Handles{nSD} = [hp_pos hp_pos hp_neg hp_neg]; %Handles of objects to animate
    animStruct.Props{nSD} = {'Vertices', 'CData', 'Vertices', 'CData'}; %Properties of objects to animate
    animStruct.Set{nSD} = {posV(:,:,nSD), posC(:,nSD), negV(:,:,nSD), negC(:,nSD)}; %Property values for to set in order to animate
    
end

%Animate figure
anim8(hf, animStruct);

%%%% TODO: could provide alternative views across subplots in anim8 figure?

%% Assess relationships between tibia & trabecular shape models

%Change directory
cd('..\ShapeModels\tibia\shape_model');

%Load matlab version
tibiaShapeModel = load('tibia.pc.mat');

%Change directory
cd('..\..\trabecular\shape_model');

%Load matlab version
trabecularShapeModel = load('trabecular.pc.mat');

%Reshape mean
tibiaShapeModel.meanPts = transpose(reshape(tibiaShapeModel.mean, [3, length(tibiaShapeModel.mean)/3]));
trabecularShapeModel.meanPts = transpose(reshape(trabecularShapeModel.mean, [3, length(trabecularShapeModel.mean)/3]));

%Normalise each row of projected weights to get them in normalised/SD form
for nPC = 1:length(tibiaShapeModel.weights)
    tibiaShapeModel.normProjectedWeights(nPC,:) = zscore(tibiaShapeModel.projectedWeights(nPC,:));
end
for nPC = 1:length(trabecularShapeModel.weights)
    trabecularShapeModel.normProjectedWeights(nPC,:) = zscore(trabecularShapeModel.projectedWeights(nPC,:));
end

% % % %Example scatter plot of relationship between correspondent PCs
% % % nPC = 2;
% % % scatter(tibiaShapeModel.normProjectedWeights(1,:), ...
% % %     trabecularShapeModel.normProjectedWeights(3,:))
% % % axis equal
% % % %Calculate corresponding correlation coefficient
% % % [r,p] = corrcoef(tibiaShapeModel.normProjectedWeights(1,:), ...
% % %     trabecularShapeModel.normProjectedWeights(3,:));


%%% PC1 to PC1 high relationship given both size factors
%%% PC2 to PC2 not the same, given they are different things
%%% PC1 from tibia to PC2 of trabecular has a reasonable relationship IF
%%% the dodgy PC2 one (ind = 9) is removed (i.e. need to fix those surfaces...)

%Is it possible to create a linear regression model from the tibia surface
%PCs to predict the individual trabecular PCs?

%Regression model for trabecular PC

%Create prediction matrix
%Use first 6 PCs as test
X = transpose(tibiaShapeModel.normProjectedWeights(1:6,:));

%Create output matrix
y = transpose(trabecularShapeModel.normProjectedWeights(3,:));

%Would removing outliers > 1 or 2 SD assist in creating a better model???
%2 doesn't make an improvement for this PC, dropping down to 1 makes it a
%better value... This isn't necessary for a 'good PC' like PC1 though - so
%maybe it's more about fixing the model? 
% X = X(find(abs(y) < 2),:);
% y = y(find(abs(y) < 2),:);

%Create linear model
mdl = fitlm(X,y)
%%%% Gives pretty good R-squared value for PC1...
%%%% Fiarly good result to for PC2 once dodgy value (ind = 9) removed...

%Plot model
plot(mdl);

%Output predicted values
yPred = predict(mdl, X);

%Calculate the root mean square error between the two
rmse = sqrt(sum((y - yPred).^2) / length(y));

%%%% End idea would be to reconstruct the trabecula from the combined
%%%% linear regression models predicted PC scores (i.e. a normalised SD
%%%% matrix of PCs) and compare this to the actual participants
%%%% surface...

%%%% The outliers for the various PCs need to be fixed up, as they cause
%%%% much poorer model fits...
%%%%
%%%% Hard to tell right now if it is the poor trabecular shape model
%%%% resulting in the secondary PCs being poorly predicted, or whether it
%%%% is just more difficult to predict these smaller variations???



%% Heat map only

%%%% Basic stuff at the moment...

%Read in created points from gias2 reconstruction
meanPoints = readmatrix('pointCloud_mean.csv');
for cc = 1:10
    
    %Load the reconstructed component points
    p3_Points = readmatrix(['pointCloud_pc',num2str(cc),'_p3.csv']);
    m3_Points = readmatrix(['pointCloud_pc',num2str(cc),'_m3.csv']);

    %%%%%%%% TODO: better looping...

    %The reconstructed points don't align with the STL points after merging,
    %but theoretically we should be able to find the matching points based on
    %the minimum distance. We calculate the distance between each individual
    %point in the mean imported point cloud and take the minimum to be the
    %matching point. This gives us the index of the row to match the data to
    %across the mean and reconstructed points so that they then match the STL.
    for pp = 1:length(meanPoints)

        %Calculate current point distances to STL nodes
        distCurrPoint = distancePoints3d(meanV(pp,:),meanPoints);

        %Get the index of the min distance points
        minInd = find(distCurrPoint == min(distCurrPoint));

        %Place the points from the reconstructed point clouds at the
        %appropriate row to match the STL
        %Note that we don't really need to reallign the mean points here as
        %they already exist from the STL. We can check it here for sanity
        %though if we want...
    % % %     meanCheck(pp,:) = meanPoints(minInd,:);
        pcPlusPoints(pp,:) = p3_Points(minInd,:);
        pcMinusPoints(pp,:) = m3_Points(minInd,:);

    end
    clear pp

% % %     %Visualise
% % %     cFigure;
% % %     gpatch(meanF,pcPlusPoints);
% % %     axisGeom;

    %Calculate distances between points for heat mapping
    pcPlusDist = distancePoints3d(meanV,pcPlusPoints);
    pcMinusDist = distancePoints3d(meanV,pcMinusPoints);

    %Convert distances to face measures
    pcPlusC = vertexToFaceMeasure(meanF,pcPlusDist);
    pcMinusC = vertexToFaceMeasure(meanF,pcMinusDist);

    %Visualise surfaces and change in points

    %Create & export figure
    
    %Lateral
    cFigure;
    gpatch(meanF,pcPlusPoints,pcPlusC,'none');
    colormap viridis; %colorbar;
    axisGeom; 
    view(90,90); %lateral
    axis off;
    camlight headlight
    export_fig(['heatmap_PC',num2str(cc),'_lateral.png'],'-m2.5')
    close()
    
    %Medial
    cFigure;
    gpatch(meanF,pcPlusPoints,pcPlusC,'none');
    colormap viridis; %colorbar;
    axisGeom; 
    view(90,-90); %medial
    axis off;
    camlight headlight
    export_fig(['heatmap_PC',num2str(cc),'_medial.png'],'-m2.5')
    close()
    
    %Posterior
    cFigure;
    gpatch(meanF,pcPlusPoints,pcPlusC,'none');
    colormap viridis; %colorbar;
    axisGeom; 
    view(90,180); %posterior
    axis off;
    camlight headlight
    export_fig(['heatmap_PC',num2str(cc),'_posterior.png'],'-m2.5')
    close()
    
    %Anterior
    cFigure;
    gpatch(meanF,pcPlusPoints,pcPlusC,'none');
    colormap viridis; %colorbar;
    axisGeom; 
    view(90,0); %anterior
    axis off;
    camlight headlight
    export_fig(['heatmap_PC',num2str(cc),'_anterior.png'],'-m2.5')
    close()

end
clear cc
