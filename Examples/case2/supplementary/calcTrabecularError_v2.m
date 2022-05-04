function [sumError] = calcTrabecularError_v2(x, shapeModel, origPoints)

    %% This script serves as the function input to the fmincon procedure
    %  where we aim to identify the PC scores that minimise the point error
    %  in the reconstruction of the a surface. The function
    %  reconstructs the surface and calculates the error
    
    %  Inputs:
    %       x - 15 x 1 array of:
    %           > 3 x 1 array of a surface scaling factor
    %           > 3 x 1 array of a translational matrix
    %           > 9 x 1 array of a rotational matrix (each interval of 3 equates to a row of the rotation matrix)
    %       shapeModel - the structure with shape model data for reconstructing surface
    %       origPoints - n x 3 array of the XYZ original points to compare to  
    
    %% Set-up
    
    %Extract the scaling factor
    scaleFactor(1:3,1) = x(1:3);
    
    %Extract the translational matrix
    transMat(1:3,1) = x(4:6);
    
    %Extract the rotation matrix
    rotMat(1,1:3) = x(7:9);
    rotMat(2,1:3) = x(10:12);
    rotMat(3,1:3) = x(13:15);
    
    %% Manipulate surface
    
    %Apply translation, rotation and scaling
    finalPoints = scaleFactor' .* origPoints * rotMat' + ...
        repmat(transMat',[length(origPoints) 1]);
        
    %Visualise original vs. reconstructed
% % %     cFigure; hold on;
% % %     gpatch(shapeModel.F, origPoints, 'gw', 'k', 0.3);
% % %     gpatch(shapeModel.F, finalPoints, 'bw', 'k', 1);
% % %     axisGeom;

    %Calculate error between reconstructed and original surface
    reconError = distancePoints3d(origPoints, finalPoints);
    
    %Convert error to colour map for visualisation
    errorColF = vertexToFaceMeasure(shapeModel.F, reconError);
    errorColV = faceToVertexMeasure(shapeModel.F, origPoints, errorColF);
    
% % %     %Visualise error
% % %     cFigure; hold on;
% % %     hpOrig = gpatch(shapeModel.F, origPoints, [200/255 200/255 200/255], 'none', 0.5);
% % %     hp = gpatch(shapeModel.F, reconstructedV, errorColV, 'none', 1);
% % %     hp.FaceColor = 'Interp'; colormap viridis
% % %     axis equal; axis tight; axis off
% % %     view(0,90);
% % %     colorbar; camlight headlight
    
    %Calculate mean error across points
    sumError = sum(reconError);
    
end