function [modelLandmarks, modelLandmarkNames, modelLandmarkInds] = createBoneBorderLandmarks_shapeModel(tibiaFibulaShapeModel, modelLandmarks, modelLandmarkNames, palpableLandmarks, nTibiaPoints, nFibulaPoints)

    %% This script serves to create the specific landmarks on the tibial and
    %  fibula borders based on the digitised points along these regions
    %
    %  Inputs:
    %       tibiaFibulaShapeModel - structure containing shape model data
    %       modelLandmarks - structure containing coordinates of model landmarks
    %       modelLandmarkNames - cell array of the model landmark names
    %       palpableLandmarks - cell array of the palpable landmark names
    %       tibiaFibulaShapeModel - structure containing tibia-fibula shape model data
    %       nTibiaPoints - number of points in shape model corresponding to tibia nodes
    %       nFibulaPoints - number of points in shape model corresponding to fibula nodes
    
    %% Create tibial landmarks
    
    %Identify landmark names for tibial border
    for landmarkNo = 1:length(modelLandmarkNames)
        if startsWith(modelLandmarkNames{landmarkNo}(1), 'T') && isstrprop(modelLandmarkNames{landmarkNo}(2),'digit')
            isBorderLandmark(landmarkNo,1) = true;
        else
            isBorderLandmark(landmarkNo,1) = false;
        end
    end
    antBorderLandmarkNames = modelLandmarkNames(isBorderLandmark);

    %Set distance points along border to create landmarks
    antBorderDists = [0.25, 0.50, 0.75];

    %Measure longitudinal axis distance between medial condyle and malleolus
    MC_MM = modelLandmarks.MC(2) - modelLandmarks.MM(2);

    %Loop through distances to calculate
    for distNo = 1:length(antBorderDists)
        %Identify Y-axis point that corresponds to current distance
        yPoint = modelLandmarks.MC(2) - (antBorderDists(distNo) * MC_MM);
        %Identify closest tibial anterior border landmarks to the Y-point
        for landmarkNo = 1:length(antBorderLandmarkNames)
            antBorderLandmarkDist(landmarkNo,1) = modelLandmarks.(antBorderLandmarkNames{landmarkNo})(2) - yPoint;
        end
        %Identify the point of positive to negative switch. Landmark will lie
        %between those points
        posNegSwitch = find(antBorderLandmarkDist(1:end-1) > 0 & antBorderLandmarkDist(2:end) < 0);
        %Create points (n = 100) on the line between these two points
        linePoints = (1-(linspace(0,1,100)')) * modelLandmarks.(antBorderLandmarkNames{posNegSwitch}) + ...
            (linspace(0,1,100)') * modelLandmarks.(antBorderLandmarkNames{posNegSwitch+1});
        %Identify the point on this line that has the closest distance to the y-point
        minDistInd = find(min(abs(yPoint - linePoints(:,2))) == abs(yPoint - linePoints(:,2)));
        %Find the closest vertex on the model surface to this to allocate the landmark
        landmarkPtInd = find(min(distancePoints3d(tibiaFibulaShapeModel.meanPoints(1:nTibiaPoints,:),linePoints(minDistInd,:))) == ...
            distancePoints3d(tibiaFibulaShapeModel.meanPoints(1:nTibiaPoints,:),linePoints(minDistInd,:)));
        %Allocate landmark to structure
        newLandmarkName = ['T',num2str(round(antBorderDists(distNo)*100))];
        modelLandmarks.(char(newLandmarkName)) = tibiaFibulaShapeModel.meanPoints(landmarkPtInd,:);
    end

    %% Create fibula landmark

    %Identify landmark names for fibula shaft
    for landmarkNo = 1:length(modelLandmarkNames)
        if startsWith(modelLandmarkNames{landmarkNo}(1), 'F') && isstrprop(modelLandmarkNames{landmarkNo}(2),'digit')
            isFibulaLandmark(landmarkNo,1) = true;
        else
            isFibulaLandmark(landmarkNo,1) = false;
        end
    end
    fibulaShaftLandmarkNames = modelLandmarkNames(isFibulaLandmark);

    %Measure longitudinal distance between fibula head and malleolus
    FH_LM = modelLandmarks.FH(2) - modelLandmarks.MM(2);

    %Calculate for 75% distance
    %Identify Y-axis point that corresponds to current distance
    yPoint = modelLandmarks.FH(2) - (0.75 * FH_LM);
    %Identify closest fibula shaft landmarks to the Y-point
    for landmarkNo = 1:length(fibulaShaftLandmarkNames)
        fibShaftLandmarkDist(landmarkNo,1) = modelLandmarks.(fibulaShaftLandmarkNames{landmarkNo})(2) - yPoint;
    end
    %Identify the point of positive to negative switch. Landmark will lie
    %between those points
    posNegSwitch = find(fibShaftLandmarkDist(1:end-1) > 0 & fibShaftLandmarkDist(2:end) < 0);
    %Create points (n = 100) on the line between these two points
    linePoints = (1-(linspace(0,1,100)')) * modelLandmarks.(fibulaShaftLandmarkNames{posNegSwitch}) + ...
        (linspace(0,1,100)') * modelLandmarks.(fibulaShaftLandmarkNames{posNegSwitch+1});
    %Identify the point on this line that has the closest distance to the y-point
    minDistInd = find(min(abs(yPoint - linePoints(:,2))) == abs(yPoint - linePoints(:,2)));
    %Find the closest vertex on the model surface to this to allocate the landmark
    landmarkPtInd = find(min(distancePoints3d(tibiaFibulaShapeModel.meanPoints(nTibiaPoints+1:end,:),linePoints(minDistInd,:))) == ...
        distancePoints3d(tibiaFibulaShapeModel.meanPoints(nTibiaPoints+1:end,:),linePoints(minDistInd,:))) + nTibiaPoints;
    %Allocate landmark to structure
    modelLandmarks.F75 = tibiaFibulaShapeModel.meanPoints(landmarkPtInd,:);

    %% Identify appropriate landmark indices
    
    %Create index labels for landmarks on the shape model. This matches up the
    %landmarks to data points in the shape model.
    %Identify closest points to those provided to identify index
    for landmarkNo = 1:length(palpableLandmarks)
        %Calculate distance between points
        ptDist = distancePoints3d(tibiaFibulaShapeModel.meanPoints, ...
            ones(length(tibiaFibulaShapeModel.meanPoints),1) * modelLandmarks.(palpableLandmarks{landmarkNo}));
        %Find index of minimum and allocate to variable
        modelLandmarkInds.(palpableLandmarks{landmarkNo}) = find(min(ptDist) == ptDist);
    end
    
    %% Review landmarks
    
% % %     %Visualise mean shape model and palpable landmarks
% % %     cFigure; hold on
% % %     gpatch(tibiaFibulaShapeModel.F,tibiaFibulaShapeModel.meanPoints,'gw','k', 0.5)
% % %     for landmarkNo = 1:length(palpableLandmarks)
% % %         plotV(modelLandmarks.(palpableLandmarks{landmarkNo}), 'r.','MarkerSize', 25);
% % %     end
% % %     axisGeom; camlight headlight    
    
end