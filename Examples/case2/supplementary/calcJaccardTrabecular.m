function [jaccardSimilarity] = calcJaccardTrabecular(actualV, predictedV, shapeModel)

    %% Function to calculate the Jaccard index - i.e. a measure of volume
    %  consistency - between the actual trabecular and predicted trabecular
    
        
    %  Inputs:
    %       actualV - n x 3 array of the XYZ points of the actual trabecular
    %       predictedV - n x 3 array of the XYZ points of the predicted trabecular
    %       shapeModel - the structure with shape model data for the trabecular    
    
    %% Calculate Jaccard similarity

    %Identify size for a blank logical mask to fill with the actual volume locations
    %Size of mask is based on maximum size of trabecular in X and Z directions
    minX = min([floor(min(actualV(:,1))),floor(min(predictedV(:,1)))]);
    maxX = max([ceil(max(actualV(:,1))),ceil(max(predictedV(:,1)))]);
    diffX = maxX - minX;
    minZ = min([floor(min(actualV(:,3))),floor(min(predictedV(:,3)))]);
    maxZ = max([ceil(max(actualV(:,3))),ceil(max(predictedV(:,3)))]);
    diffZ = maxZ - minZ;

    %Create look-up variables of X and Z coordinate steps to look trhough and
    %apply to mask
    xCoordSteps = linspace(minX, maxX, diffX+1);
    zCoordSteps = linspace(minZ, maxZ, diffZ+1);

    %Create the set of y-coordinates to look through based on min and max coordinates
    minY = min([floor(min(actualV(:,2))),floor(min(predictedV(:,2)))]);
    maxY = max([ceil(max(actualV(:,2))),ceil(max(predictedV(:,2)))]);
    diffY = maxY - minY;
    yCoordSteps = linspace(minY, maxY, diffY+1);
    yCoordSteps = yCoordSteps(1:2:end);

    %Set starting values for TP, FP and FN
    TP = 0; FP = 0; FN = 0;
    
    %Create wait bar for calculating Jaccard similarity
    wbar = waitbar(0, 'Calculating Jaccard similarity...');
    
    %Loop through y-coordinate steps and calculate Jaccard similarity
    for yInd = 1:length(yCoordSteps)

        %Set y-level to search on current iteration
        yLevel = yCoordSteps(yInd);

        %Create set of points to check based on current coordinate steps
        ptInd = 1;
        for row = 1:length(xCoordSteps)
            for col = 1:length(zCoordSteps)
                chkPts(ptInd,:) = [xCoordSteps(row), yLevel, zCoordSteps(col)];
                ptInd = ptInd+1;
            end
        end

        %Check whether points are within original and reconsructed surfaces
        origInside = intriangulation(actualV, shapeModel.F, chkPts);
        recInside = intriangulation(predictedV, shapeModel.F, chkPts);

        %Reshape logical values to fit mask shape
        maskOrig = reshape(origInside', [], length(xCoordSteps))';
        maskRec = reshape(recInside', [], length(xCoordSteps))';

        %Manually calculate jaccard similarity
        TP = TP + sum(maskOrig & maskRec, 'all');
        FP = FP + sum(maskOrig & ~maskRec, 'all');
        FN = FN + sum(~maskOrig & maskRec, 'all');

    % % %     %Compare masks
    % % %     figure;
    % % %     imshowpair(maskOrig, maskRec);
    
        %Update waitbar
        waitbar(yInd / length(yCoordSteps), wbar);

    end

    %Calculate the jaccard similarity of the total slices
    jaccardSimilarity = TP / (TP + FP + FN);
    
    %Close waitbar
    close(wbar);
    
end