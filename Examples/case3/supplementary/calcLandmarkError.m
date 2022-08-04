function [sumError] = calcLandmarkError(pcScores, shapeModel, sampleLandmarks, shapeModelLandmarkInds)

    %% This script serves as the function input to the fmincon procedure
    %  where we aim to identify the PC scores that minimise the point error
    %  in the reconstruction of the a surface. The function
    %  reconstructs the surface and calculates the error
    
    %  Inputs:
    %       pcScores - n x 1 array of PC scores to use in reconstruction
    %       shapeModel - the structure with shape model data for reconstructing surface
    %       landmarks - n x 3 array of the XYZ original landmarks to compare shape model to  
    
    %% Set-up
    
    %Calculate the number of PCs to use based on the score input
    nPC = length(pcScores);
    
    %% Reconstruct surface
    
    %Use the provided scores and shape model loadings to reconstruct surface
    %Reconstruct the tibia points for new surface
    reconstructedPoints = pcScores' * ...
        shapeModel.loadings(:,1:nPC)' + ...
        shapeModel.mean;

    %Reshape points to visualise
    reconstructedV = transpose(reshape(reconstructedPoints,...
        [3, length(reconstructedPoints)/3]));
    
    %Extract shape model landmark locations
    landmarkNames = fieldnames(shapeModelLandmarkInds);
    for landmarkNo = 1:length(landmarkNames)
        reconstructedLandmarks.(landmarkNames{landmarkNo}) = ...
            reconstructedV(shapeModelLandmarkInds.(landmarkNames{landmarkNo}),:);        
    end
        
% % %     %Visualise original vs. reconstructed and landmarks
% % %     cFigure; hold on;
% % %     %Reconstructed shape
% % %     gpatch(shapeModel.F, reconstructedV, 'bw', 'none', 0.3);
% % %     %Reconstructed & sample landmarks
% % %     for landmarkNo = 1:length(landmarkNames)
% % %         plotV(reconstructedLandmarks.(landmarkNames{landmarkNo}), 'r.', 'MarkerSize', 25)
% % %         plotV(sampleLandmarks.(landmarkNames{landmarkNo}), 'g.', 'MarkerSize', 25)
% % %     end
% % %     axisGeom; camlight headlight
    
    %% Realign shape model reconstruction to original landmarks
    
    %Extract landmarks into array
    for landmarkNo = 1:length(landmarkNames)
        reconLandmarks(landmarkNo,:) = reconstructedLandmarks.(landmarkNames{landmarkNo});
        origLandmarks(landmarkNo,:) = sampleLandmarks.(landmarkNames{landmarkNo});
    end
    
    %Rigidly align the reconstructed landmarks to the original
    [~, Z] = procrustes(origLandmarks, reconLandmarks,...
        'Scaling', false, 'Reflection', false);
    
    %Extract out the transformed landmarks to a named structure
    for landmarkNo = 1:length(landmarkNames)
        reconstructedAlignedLandmarks.(landmarkNames{landmarkNo}) = Z(landmarkNo,:);
    end
    
    %% Calculate landmark error
    
    %Calculate error between sample landmarks and shape model landmarks
    for landmarkNo = 1:length(landmarkNames)
        reconError(landmarkNo) = distancePoints3d(sampleLandmarks.(landmarkNames{landmarkNo}), ...
            reconstructedAlignedLandmarks.(landmarkNames{landmarkNo}));
    end
    
    %Calculate sum of error across points
    sumError = sum(reconError);
    
end