function [sampleLandmarks, sampleLandmarkNames] = createBoneBorderLandmarks_sample(sampleLandmarks, sampleLandmarkNames, palpableLandmarks, tibiaV, fibulaV)

    %% This script serves to create the specific landmarks on the tibial and
    %  fibula borders based on the digitised points along these regions
    %
    %  Inputs:
    %       modelLandmarks - structure containing coordinates of model landmarks
    %       modelLandmarkNames - cell array of the model landmark names
    %       palpableLandmarks - cell array of the palpable landmark names
    %       tibiaV - structure containing tibia node data
    %       fibulaV - structure containing fibula node data
    
    %% Create tibial landmarks
    
    %Identify landmark names for tibial border
    for landmarkNo = 1:length(sampleLandmarkNames)
        if startsWith(sampleLandmarkNames{landmarkNo}(1), 'T') && isstrprop(sampleLandmarkNames{landmarkNo}(2),'digit')
            isBorderLandmark(landmarkNo,1) = true;
        else
            isBorderLandmark(landmarkNo,1) = false;
        end
    end
    antBorderLandmarkNames = sampleLandmarkNames(isBorderLandmark);
    
    %Set distance points along border to create landmarks
    antBorderDists = [0.25, 0.50, 0.75];

    %Measure longitudinal axis distance between medial condyle and malleolus
    MC_MM = sampleLandmarks.MC(2) - sampleLandmarks.MM(2);

    %Loop through distances to calculate
    for distNo = 1:length(antBorderDists)
        %Identify Y-axis point that corresponds to current distance
        yPoint = sampleLandmarks.MC(2) - (antBorderDists(distNo) * MC_MM);
        %Identify closest tibial anterior border landmarks to the Y-point
        for landmarkNo = 1:length(antBorderLandmarkNames)
            antBorderLandmarkDist(landmarkNo,1) = sampleLandmarks.(antBorderLandmarkNames{landmarkNo})(2) - yPoint;
        end
        %Identify the point of positive to negative switch. Landmark will lie
        %between those points
        posNegSwitch = find(antBorderLandmarkDist(1:end-1) > 0 & antBorderLandmarkDist(2:end) < 0);
        %Create points (n = 100) on the line between these two points
        linePoints = (1-(linspace(0,1,100)')) * sampleLandmarks.(antBorderLandmarkNames{posNegSwitch}) + ...
            (linspace(0,1,100)') * sampleLandmarks.(antBorderLandmarkNames{posNegSwitch+1});
        %Identify the point on this line that has the closest distance to the y-point
        minDistInd = find(min(abs(yPoint - linePoints(:,2))) == abs(yPoint - linePoints(:,2)));
        %Find the closest vertex on the model surface to this to allocate the landmark
        landmarkPtInd = find(min(distancePoints3d(tibiaV,linePoints(minDistInd,:))) == ...
            distancePoints3d(tibiaV,linePoints(minDistInd,:)));
        %Allocate landmark to structure
        newLandmarkName = ['T',num2str(round(antBorderDists(distNo)*100))];
        sampleLandmarks.(char(newLandmarkName)) = tibiaV(landmarkPtInd,:);
    end

    %% Create fibula landmark

    %Identify landmark names for fibula shaft
    for landmarkNo = 1:length(sampleLandmarkNames)
        if startsWith(sampleLandmarkNames{landmarkNo}(1), 'F') && isstrprop(sampleLandmarkNames{landmarkNo}(2),'digit')
            isFibulaLandmark(landmarkNo,1) = true;
        else
            isFibulaLandmark(landmarkNo,1) = false;
        end
    end
    fibulaShaftLandmarkNames = sampleLandmarkNames(isFibulaLandmark);

    %Measure longitudinal distance between fibula head and malleolus
    FH_LM = sampleLandmarks.FH(2) - sampleLandmarks.MM(2);

    %Calculate for 75% distance
    %Identify Y-axis point that corresponds to current distance
    yPoint = sampleLandmarks.FH(2) - (0.75 * FH_LM);
    %Identify closest fibula shaft landmarks to the Y-point
    for landmarkNo = 1:length(fibulaShaftLandmarkNames)
        fibShaftLandmarkDist(landmarkNo,1) = sampleLandmarks.(fibulaShaftLandmarkNames{landmarkNo})(2) - yPoint;
    end
    %Identify the point of positive to negative switch. Landmark will lie
    %between those points
    posNegSwitch = find(fibShaftLandmarkDist(1:end-1) > 0 & fibShaftLandmarkDist(2:end) < 0);
    %Create points (n = 100) on the line between these two points
    linePoints = (1-(linspace(0,1,100)')) * sampleLandmarks.(fibulaShaftLandmarkNames{posNegSwitch}) + ...
        (linspace(0,1,100)') * sampleLandmarks.(fibulaShaftLandmarkNames{posNegSwitch+1});
    %Identify the point on this line that has the closest distance to the y-point
    minDistInd = find(min(abs(yPoint - linePoints(:,2))) == abs(yPoint - linePoints(:,2)));
    %Find the closest vertex on the model surface to this to allocate the landmark
    landmarkPtInd = find(min(distancePoints3d(fibulaV,linePoints(minDistInd,:))) == ...
        distancePoints3d(fibulaV,linePoints(minDistInd,:)));
    %Allocate landmark to structure
    sampleLandmarks.F75 = fibulaV(landmarkPtInd,:);
    
end