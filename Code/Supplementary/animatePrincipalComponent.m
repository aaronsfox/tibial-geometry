function animatePrincipalComponent(shapeModel, principalComponent, sdRange, exportAnimation)

    %% This function creates an animation of the shape change caused by
    %  manipulating a principal component score
    %
    %  Inputs:
    %
    %   shapeModel - the data structure containing shape model details
    %   principalComponent - the PC number to manipulate
    %   sdRange - value representing +/- SD range to visualise (default = 3)
    %   exportAnimation - option to export the created animations (default = true)

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
    if nargin < 4
        exportAnimation = true
    end
    
    %Set array to rotate views of surfaces later in reconstruction process
    surfaceRot = [-90, 0, 90, 180];

    %Set list to label reconstruction subplots with
    viewLabel = [{'Anterior'}, {'Lateral'}, {'Posterior'}, {'Medial'}];
    
    %% Reconstruct surfaces across steps of the SD range
    
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
    
    %% Create the plus SD figure
    
    %Create the basic view to store graphics and initiate animation
    hfPos = cFigure; hold on;
    
    %Loop through four views to create subplot
    for viewNo = 1:4    
        %Create subplot for current view
        subplot(1,4,viewNo);
        %Add patch data
        posFigHandles.(viewLabel{viewNo}).hp_pos = gpatch(shapeModel.F, meanPts, posC(:,1), 'none', 1);
        %Interpolate colouring for smoothness
        posFigHandles.(viewLabel{viewNo}).hp_pos.FaceColor = 'Interp'; colormap viridis
        %Set axes parameters
        axis equal; axis tight;
        view(0,90); rotate(posFigHandles.(viewLabel{viewNo}).hp_pos,[0 1 0], surfaceRot(viewNo));
        %Set axis parameters
        camlight headlight; axis off
        %Add title
        title([viewLabel{viewNo},' View'], 'FontSize', 12);
        %Set axis limits based on data values
        for nSD = 1:length(calcSD)
            posLim(:,:,nSD) = transformPoint3d(posV(:,:,nSD),createRotationOy(deg2rad(surfaceRot(viewNo))));
            negLim(:,:,nSD) = transformPoint3d(negV(:,:,nSD),createRotationOy(deg2rad(surfaceRot(viewNo))));
        end
        axis(axisLim([posLim;negLim]));
    end
    
    %Set up animation
    animStruct.Time = calcSD;
    %Loop through SD points
    for nSD = 1:length(calcSD)

        %Set entries into animation structure
        %Handles of objects to animate
        animStruct.Handles{nSD} = [posFigHandles.(viewLabel{1}).hp_pos posFigHandles.(viewLabel{1}).hp_pos ...
            posFigHandles.(viewLabel{2}).hp_pos posFigHandles.(viewLabel{2}).hp_pos ...
            posFigHandles.(viewLabel{3}).hp_pos posFigHandles.(viewLabel{3}).hp_pos ...
            posFigHandles.(viewLabel{4}).hp_pos posFigHandles.(viewLabel{4}).hp_pos];
        %Properties of objects to animate
        animStruct.Props{nSD} = {'Vertices', 'CData', ...
            'Vertices', 'CData', ...
            'Vertices', 'CData', ...
            'Vertices', 'CData'};
        %Property values for to set in order to animate
        animStruct.Set{nSD} = {transformPoint3d(posV(:,:,nSD),createRotationOy(deg2rad(surfaceRot(1)))), posC(:,nSD), ...
            transformPoint3d(posV(:,:,nSD),createRotationOy(deg2rad(surfaceRot(2)))) posC(:,nSD), ...
            transformPoint3d(posV(:,:,nSD),createRotationOy(deg2rad(surfaceRot(3)))), posC(:,nSD), ...
            transformPoint3d(posV(:,:,nSD),createRotationOy(deg2rad(surfaceRot(4)))), posC(:,nSD)}; 

    end
    
    %Set labelling option
    animStruct.labelOption = 'SD';
    
    %Add overall figure title
    sgtitle(['Surface Shape Change from PC',num2str(principalComponent),' with +SD'],...
        'FontSize', 16, 'FontWeight', 'bold')
    
    %Animate figure
    anim8_labelOption(hfPos, animStruct);
    
    %Export if desired
    if exportAnimation
        
        %Export animation
        exportGifAnim8();
        
        %Unpack from export folder and rename
        %Get name of file in export folder
        fDir = dir('efw');
        %Move file
        movefile(['efw\',fDir(end).name], ...
            ['PC',num2str(principalComponent),'_plus-',num2str(sdRange),'SD_animation.gif']);
        %Remove export directory
        rmdir('efw');
        
    end
    
    %Close figure
    close(hfPos);
    
    %% Create the negative SD figure
    
    %Create the basic view to store graphics and initiate animation
    hfNeg = cFigure; hold on;
    
    %Loop through four views to create subplot
    for viewNo = 1:4    
        %Create subplot for current view
        subplot(1,4,viewNo);
        %Add patch data
        negFigHandles.(viewLabel{viewNo}).hp_neg = gpatch(shapeModel.F, meanPts, posC(:,1), 'none', 1);
        %Interpolate colouring for smoothness
        negFigHandles.(viewLabel{viewNo}).hp_neg.FaceColor = 'Interp'; colormap viridis
        %Set axes parameters
        axis equal; axis tight;
        view(0,90); rotate(negFigHandles.(viewLabel{viewNo}).hp_neg,[0 1 0], surfaceRot(viewNo));
        %Set axis parameters
        camlight headlight; axis off
        %Add title
        title([viewLabel{viewNo},' View'], 'FontSize', 12);
        %Set axis limits based on data values
        for nSD = 1:length(calcSD)
            posLim(:,:,nSD) = transformPoint3d(posV(:,:,nSD),createRotationOy(deg2rad(surfaceRot(viewNo))));
            negLim(:,:,nSD) = transformPoint3d(negV(:,:,nSD),createRotationOy(deg2rad(surfaceRot(viewNo))));
        end
        axis(axisLim([posLim;negLim]));
    end
    
    %Set up animation
    animStruct.Time = calcSD;
    %Loop through SD points
    for nSD = 1:length(calcSD)

        %Set entries into animation structure
        %Handles of objects to animate
        animStruct.Handles{nSD} = [negFigHandles.(viewLabel{1}).hp_neg negFigHandles.(viewLabel{1}).hp_neg ...
            negFigHandles.(viewLabel{2}).hp_neg negFigHandles.(viewLabel{2}).hp_neg ...
            negFigHandles.(viewLabel{3}).hp_neg negFigHandles.(viewLabel{3}).hp_neg ...
            negFigHandles.(viewLabel{4}).hp_neg negFigHandles.(viewLabel{4}).hp_neg];
        %Properties of objects to animate
        animStruct.Props{nSD} = {'Vertices', 'CData', ...
            'Vertices', 'CData', ...
            'Vertices', 'CData', ...
            'Vertices', 'CData'};
        %Property values for to set in order to animate
        animStruct.Set{nSD} = {transformPoint3d(negV(:,:,nSD),createRotationOy(deg2rad(surfaceRot(1)))), negC(:,nSD), ...
            transformPoint3d(negV(:,:,nSD),createRotationOy(deg2rad(surfaceRot(2)))) negC(:,nSD), ...
            transformPoint3d(negV(:,:,nSD),createRotationOy(deg2rad(surfaceRot(3)))), negC(:,nSD), ...
            transformPoint3d(negV(:,:,nSD),createRotationOy(deg2rad(surfaceRot(4)))), negC(:,nSD)}; 

    end
    
    %Set labelling option
    animStruct.labelOption = 'SD';
    
    %Add overall figure title
    sgtitle(['Surface Shape Change from PC',num2str(principalComponent),' with -SD'],...
        'FontSize', 16, 'FontWeight', 'bold')
    
    %Animate figure
    anim8_labelOption(hfNeg, animStruct);
    
    %Export if desired
    if exportAnimation
        
        %Export animation
        exportGifAnim8();
        
        %Unpack from export folder and rename
        %Get name of file in export folder
        fDir = dir('efw');
        %Move file
        movefile(['efw\',fDir(end).name], ...
            ['PC',num2str(principalComponent),'_minus-',num2str(sdRange),'SD_animation.gif']);
        %Remove export directory
        rmdir('efw');
        
    end
    
    %Close figure
    close(hfNeg);
    
end