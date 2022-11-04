function animatePrincipalComponentTwoParts(part1_nPoints, part2_nPoints, shapeModel, principalComponent, sdRange, exportAnimation)

    %% This function creates an animation of the shape change caused by
    %  manipulating a principal component score. This version allows a
    %  visualisation of two parts. It's assumed that part 1 is on the
    %  outside and part 2 on the inside. It's also assumed the points are
    %  arranged in order of part 1 then part 2. There is also an assumption
    %  that the shape model has two separate faces variables.
    %
    %  Inputs:
    %
    %   part1_nPoints - number of 3D points in the first part
    %   part2_nPoints - number of 3D points in the second part
    %   shapeModel - the data structure containing shape model details
    %   principalComponent - the PC number to manipulate
    %   sdRange - value representing +/- SD range to visualise (default = 3)
    %   exportAnimation - option to export the created animations (default = true)

    %Check inputs
    if nargin < 2
        %Error out as no shape model data provided
        error('The number of 3D points in each part is required.');
    end    
    if nargin < 3
        %Error out as no shape model data provided
        error('Shape model data structure required.');
    end
    if nargin < 4
        %Request input of which principal component to examine
        principalComponent = input('Enter principal component number to visualise: ');
    end
    if nargin < 5
        %Default to +/- 3SD's
        sdRange = 3;
    end
    if nargin < 6
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
        %Part 1
        pcPlusC1 = vertexToFaceMeasure(shapeModel.F1, pcPlusDist(1:part1_nPoints));
        pcMinusC1 = vertexToFaceMeasure(shapeModel.F1, pcMinusDist(1:part1_nPoints));
        %Part 2
        pcPlusC2 = vertexToFaceMeasure(shapeModel.F2, pcPlusDist(part1_nPoints+1:end));
        pcMinusC2 = vertexToFaceMeasure(shapeModel.F2, pcMinusDist(part1_nPoints+1:end));
        
        %Convert face colouring to vertices for colour interpolation
        %Part 1
        pcPlusCV1 = faceToVertexMeasure(shapeModel.F1, posV(1:part1_nPoints,:,nSD), pcPlusC1);
        pcMinusCV1 = faceToVertexMeasure(shapeModel.F1, negV(1:part1_nPoints,:,nSD), pcMinusC1);
        %Part 2
        pcPlusCV2 = faceToVertexMeasure(shapeModel.F2, posV(part1_nPoints+1:end,:,nSD), pcPlusC2);
        pcMinusCV2 = faceToVertexMeasure(shapeModel.F2, negV(part1_nPoints+1:end,:,nSD), pcMinusC2);

        %Set into structure
        %Part 1
        posC1(:,nSD) = pcPlusCV1;
        negC1(:,nSD) = pcMinusCV1; 
        %Part 2
        posC2(:,nSD) = pcPlusCV2;
        negC2(:,nSD) = pcMinusCV2; 

    end
    
    %% Create the plus SD figure
    
    %Create the basic view to store graphics and initiate animation
    hfPos = cFigure; hold on;
    
    %Loop through four views to create subplot
    for viewNo = 1:4    
        %Create subplot for current view
        subplot(1,4,viewNo);
        %Add patch data
        posFigHandles.(viewLabel{viewNo}).hp_pos1 = gpatch(shapeModel.F1, meanPts(1:part1_nPoints,:), posC1(:,1), 'none', 0.3); %part 1 see through
        posFigHandles.(viewLabel{viewNo}).hp_pos2 = gpatch(shapeModel.F2, meanPts(part1_nPoints+1:end,:), posC2(:,1), 'none', 1); %part 2 full opacity
        %Interpolate colouring for smoothness
        posFigHandles.(viewLabel{viewNo}).hp_pos1.FaceColor = 'Interp'; colormap viridis
        posFigHandles.(viewLabel{viewNo}).hp_pos2.FaceColor = 'Interp'; colormap viridis
        %Set axes parameters
        axis equal; axis tight; view(0,90);
        rotate(posFigHandles.(viewLabel{viewNo}).hp_pos1,[0 1 0], surfaceRot(viewNo));
        rotate(posFigHandles.(viewLabel{viewNo}).hp_pos2,[0 1 0], surfaceRot(viewNo));
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
        animStruct.Handles{nSD} = [posFigHandles.(viewLabel{1}).hp_pos1 posFigHandles.(viewLabel{1}).hp_pos1 ...
            posFigHandles.(viewLabel{2}).hp_pos1 posFigHandles.(viewLabel{2}).hp_pos1 ...
            posFigHandles.(viewLabel{3}).hp_pos1 posFigHandles.(viewLabel{3}).hp_pos1 ...
            posFigHandles.(viewLabel{4}).hp_pos1 posFigHandles.(viewLabel{4}).hp_pos1 ...
            posFigHandles.(viewLabel{1}).hp_pos2 posFigHandles.(viewLabel{1}).hp_pos2 ...
            posFigHandles.(viewLabel{2}).hp_pos2 posFigHandles.(viewLabel{2}).hp_pos2 ...
            posFigHandles.(viewLabel{3}).hp_pos2 posFigHandles.(viewLabel{3}).hp_pos2 ...
            posFigHandles.(viewLabel{4}).hp_pos2 posFigHandles.(viewLabel{4}).hp_pos2];
        %Properties of objects to animate
        animStruct.Props{nSD} = {'Vertices', 'CData', ...
            'Vertices', 'CData', ...
            'Vertices', 'CData', ...
            'Vertices', 'CData', ...
            'Vertices', 'CData', ...
            'Vertices', 'CData', ...
            'Vertices', 'CData', ...
            'Vertices', 'CData'};
        %Property values for to set in order to animate
        animStruct.Set{nSD} = {transformPoint3d(posV(1:part1_nPoints,:,nSD),createRotationOy(deg2rad(surfaceRot(1)))), posC1(:,nSD), ...
            transformPoint3d(posV(1:part1_nPoints,:,nSD),createRotationOy(deg2rad(surfaceRot(2)))) posC1(:,nSD), ...
            transformPoint3d(posV(1:part1_nPoints,:,nSD),createRotationOy(deg2rad(surfaceRot(3)))), posC1(:,nSD), ...
            transformPoint3d(posV(1:part1_nPoints,:,nSD),createRotationOy(deg2rad(surfaceRot(4)))), posC1(:,nSD), ...
            transformPoint3d(posV(part1_nPoints+1:end,:,nSD),createRotationOy(deg2rad(surfaceRot(1)))), posC2(:,nSD), ...
            transformPoint3d(posV(part1_nPoints+1:end,:,nSD),createRotationOy(deg2rad(surfaceRot(2)))) posC2(:,nSD), ...
            transformPoint3d(posV(part1_nPoints+1:end,:,nSD),createRotationOy(deg2rad(surfaceRot(3)))), posC2(:,nSD), ...
            transformPoint3d(posV(part1_nPoints+1:end,:,nSD),createRotationOy(deg2rad(surfaceRot(4)))), posC2(:,nSD)}; 

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
        negFigHandles.(viewLabel{viewNo}).hp_neg1 = gpatch(shapeModel.F1, meanPts(1:part1_nPoints,:), posC1(:,1), 'none', 0.3); %part 1 see through
        negFigHandles.(viewLabel{viewNo}).hp_neg2 = gpatch(shapeModel.F2, meanPts(part1_nPoints+1:end,:), posC2(:,1), 'none', 1); %part 2 full opacity
        %Interpolate colouring for smoothness
        negFigHandles.(viewLabel{viewNo}).hp_neg1.FaceColor = 'Interp'; colormap viridis
        negFigHandles.(viewLabel{viewNo}).hp_neg2.FaceColor = 'Interp'; colormap viridis
        %Set axes parameters
        axis equal; axis tight; view(0,90);
        rotate(negFigHandles.(viewLabel{viewNo}).hp_neg1,[0 1 0], surfaceRot(viewNo));
        rotate(negFigHandles.(viewLabel{viewNo}).hp_neg2,[0 1 0], surfaceRot(viewNo));
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
        animStruct.Handles{nSD} = [negFigHandles.(viewLabel{1}).hp_neg1 negFigHandles.(viewLabel{1}).hp_neg1 ...
            negFigHandles.(viewLabel{2}).hp_neg1 negFigHandles.(viewLabel{2}).hp_neg1 ...
            negFigHandles.(viewLabel{3}).hp_neg1 negFigHandles.(viewLabel{3}).hp_neg1 ...
            negFigHandles.(viewLabel{4}).hp_neg1 negFigHandles.(viewLabel{4}).hp_neg1 ...
            negFigHandles.(viewLabel{1}).hp_neg2 negFigHandles.(viewLabel{1}).hp_neg2 ...
            negFigHandles.(viewLabel{2}).hp_neg2 negFigHandles.(viewLabel{2}).hp_neg2 ...
            negFigHandles.(viewLabel{3}).hp_neg2 negFigHandles.(viewLabel{3}).hp_neg2 ...
            negFigHandles.(viewLabel{4}).hp_neg2 negFigHandles.(viewLabel{4}).hp_neg2];
        %Properties of objects to animate
        animStruct.Props{nSD} = {'Vertices', 'CData', ...
            'Vertices', 'CData', ...
            'Vertices', 'CData', ...
            'Vertices', 'CData', ...
            'Vertices', 'CData', ...
            'Vertices', 'CData', ...
            'Vertices', 'CData', ...
            'Vertices', 'CData'};
        %Property values for to set in order to animate
        animStruct.Set{nSD} = {transformPoint3d(negV(1:part1_nPoints,:,nSD),createRotationOy(deg2rad(surfaceRot(1)))), negC1(:,nSD), ...
            transformPoint3d(negV(1:part1_nPoints,:,nSD),createRotationOy(deg2rad(surfaceRot(2)))) negC1(:,nSD), ...
            transformPoint3d(negV(1:part1_nPoints,:,nSD),createRotationOy(deg2rad(surfaceRot(3)))), negC1(:,nSD), ...
            transformPoint3d(negV(1:part1_nPoints,:,nSD),createRotationOy(deg2rad(surfaceRot(4)))), negC1(:,nSD), ...
            transformPoint3d(negV(part1_nPoints+1:end,:,nSD),createRotationOy(deg2rad(surfaceRot(1)))), negC2(:,nSD), ...
            transformPoint3d(negV(part1_nPoints+1:end,:,nSD),createRotationOy(deg2rad(surfaceRot(2)))) negC2(:,nSD), ...
            transformPoint3d(negV(part1_nPoints+1:end,:,nSD),createRotationOy(deg2rad(surfaceRot(3)))), negC2(:,nSD), ...
            transformPoint3d(negV(part1_nPoints+1:end,:,nSD),createRotationOy(deg2rad(surfaceRot(4)))), negC2(:,nSD)}; 

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