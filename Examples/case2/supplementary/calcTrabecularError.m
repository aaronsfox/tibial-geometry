function [sumError] = calcTrabecularError(x, shapeModel, origPoints)

    %% This script serves as the function input to the fmincon procedure
    %  where we aim to identify the PC scores that minimise the point error
    %  in the reconstruction of the a surface. The function
    %  reconstructs the surface and calculates the error
    
    %  Inputs:
    %       x - n x 1 array of:
    %           > n x retainPCs array of PC scores to use in reconstruction
    %           > 3 x 1 array of a surface scaling factor
    %           > 3 x 1 array of a translational matrix
    %           > 9 x 1 array of a rotational matrix (each interval of 3 equates to a row of the rotation matrix)
    %       shapeModel - the structure with shape model data for reconstructing surface
    %       origPoints - n x 3 array of the XYZ original points to compare to  
    
    %% Set-up
    
    %Calculate the number of PCs to use based on the shape model input
    nPC = shapeModel.retainPCs;
    
    %Extract the shape model scores from the function inpiut
    pcScores(1:nPC,1) = x(1:nPC);
    
    %Extract the scaling factor
    scaleFactor(1:3,1) = x(nPC+1:nPC+3);
    
    %Extract the translational matrix
    transMat(1:3,1) = x(nPC+4:nPC+6);
    
    %Extract the rotation matrix
    rotMat(1,1:3) = x(nPC+7:nPC+9);
    rotMat(2,1:3) = x(nPC+10:nPC+12);
    rotMat(3,1:3) = x(nPC+13:nPC+15);
    
    %% Reconstruct surface
    
    %Use the provided scores and shape model loadings to reconstruct surface
    %Reconstruct the tibia points for new surface
    reconstructedPoints = pcScores' * ...
        shapeModel.loadings(:,1:nPC)' + ...
        shapeModel.mean;

    %Reshape points
    reconstructedV = transpose(reshape(reconstructedPoints,...
        [3, length(reconstructedPoints)/3]));
    
    %Apply translation, rotation and scaling
    finalPoints = scaleFactor' .* reconstructedV * rotMat' + ...
        repmat(transMat',[length(reconstructedV) 1]);
        
    %Visualise original vs. reconstructed
% % %     cFigure; hold on;
% % %     gpatch(shapeModel.F, origPoints, 'gw', 'k', 0.3);
% % %     gpatch(shapeModel.F, reconstructedV, 'rw', 'k', 1);
% % %     gpatch(shapeModel.F, finalPoints, 'bw', 'k', 1);
% % %     axisGeom;

    %Calculate error between reconstructed and original surface
    reconError = distancePoints3d(origPoints, finalPoints);
    
% % %     %Convert error to colour map for visualisation
% % %     errorColF = vertexToFaceMeasure(shapeModel.F, reconError);
% % %     errorColV = faceToVertexMeasure(shapeModel.F, origPoints, errorColF);
    
% % %     %Visualise error
% % %     cFigure; hold on;
% % %     hpOrig = gpatch(shapeModel.F, origPoints, [200/255 200/255 200/255], 'none', 0.5);
% % %     hp = gpatch(shapeModel.F, reconstructedV, errorColV, 'none', 1);
% % %     hp.FaceColor = 'Interp'; colormap viridis
% % %     axis equal; axis tight; axis off
% % %     view(0,90);
% % %     colorbar; camlight headlight
    
    %Calculate mean error across points
    sumError = mean(reconError);
    
end