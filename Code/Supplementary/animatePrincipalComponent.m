function animatePrincipalComponent(shapeModel, principalComponent, sdRange)

    %% This function creates an animation of the shape change caused by
    %  manipulating a principal component score
    %
    %  Inputs:
    %
    %   shapeModel - the data structure containing shape model details
    %   principalComponent - the PC number to manipulate
    %   sdRange - value representing +/- SD range to visualise (default = 3)

    %%%% TODO: use shift to vertex colouring and interpolate in animations...

    %Check inputs
    if nargin < 1
        %Error out as no shape model data provided
        error('Shape model data structure required.');
    end
    if nargin < 2
        %Request input of which principal component to examine
        principalComponent = input('Enter principal component number to visualise: ');
    end
    if nargin < 3
        %Default to +/- 3SD's
        sdRange = 3;
    end

    %Set number of PCs
    reconPCs = length(shapeModel.varExplained);

    %Reshape mean
    meanPts = reshape(shapeModel.mean, [3, length(shapeModel.mean)/3])';

    %Calculate mean and standard deviation values for scores
    sdScore = std(shapeModel.score);

    %Create a variable to work through 0.1 increments of SD range
    calcSD = linspace(0.1, sdRange, diff([0,sdRange])*10);

    %Calculate the points at each increment
    for nSD = 1:length(calcSD)

        %Reconstruct the points at the current SD interval

        %Set the simulated scores array with the standard deviation added/subtracted
        %Plus SD
        simScorePlus = zeros(1,reconPCs);
        simScorePlus(principalComponent) = simScorePlus(principalComponent) + (sdScore(principalComponent) * calcSD(nSD));
        %Minus SD
        simScoreMinus = zeros(1,reconPCs);
        simScoreMinus(principalComponent) = simScoreMinus(principalComponent) + (sdScore(principalComponent) * -calcSD(nSD));

        %Reconstruct the points, add the mean and reshape to 3D
        plusPts = reshape((simScorePlus * shapeModel.loadings(:,1:reconPCs)') + shapeModel.mean, ...
            [3, length(shapeModel.mean)/3])';
        minusPts = reshape((simScoreMinus * shapeModel.loadings(:,1:reconPCs)') + shapeModel.mean, ...
            [3, length(shapeModel.mean)/3])';

        %Set into structure
        posV(:,:,nSD) = plusPts;
        negV(:,:,nSD) = minusPts;

        %Calculate distances between points for heat mapping
        pcPlusDist = distancePoints3d(meanPts, plusPts);
        pcMinusDist = distancePoints3d(meanPts, minusPts);

        %Convert distances to face measures
        pcPlusC = vertexToFaceMeasure(shapeModel.F, pcPlusDist);
        pcMinusC = vertexToFaceMeasure(shapeModel.F, pcMinusDist);
        
        %Convert face colouring to vertices for colour interpolation
        pcPlusCV = faceToVertexMeasure(shapeModel.F, posV(:,:,nSD), pcPlusC);
        pcMinusCV = faceToVertexMeasure(shapeModel.F, negV(:,:,nSD), pcMinusC);

        %Set into structure
        posC(:,nSD) = pcPlusCV;
        negC(:,nSD) = pcMinusCV; 

    end

    %Create the basic view to store graphics and initiate animation
    hf = cFigure;
    
    %Positive SD
    %Create subplot
    subplot(1,2,1);
    %Add title
    title(['Plus SD for PC',num2str(principalComponent)]);
    %Add patch data
    hp_pos = gpatch(shapeModel.F, meanPts, posC(:,1), 'none', 1);
    %Interpolate colouring for smoothness
    hp_pos.FaceColor = 'Interp'; colormap viridis
    %Set axes parameters
    axis equal; axis tight; view(90,0);
    camlight headlight; axis off
    %Set axis limits based on data values
    axis(axisLim([posV;negV]));
    
    %Negative SD
    %Create subplot
    subplot(1,2,2);
    %Add title
    title(['Minus SD for PC',num2str(principalComponent)]);
    %Add patch data
    hp_neg = gpatch(shapeModel.F(:,:,1), meanPts, negC(:,1), 'none', 1);
    %Interpolate colouring for smoothness
    hp_neg.FaceColor = 'Interp'; colormap viridis
    %Set axes parameters
    axis equal; axis tight; view(90,0); 
    camlight headlight; axis off
    %Set axis limits based on data values
    axis(axisLim([posV;negV]));

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
    
end