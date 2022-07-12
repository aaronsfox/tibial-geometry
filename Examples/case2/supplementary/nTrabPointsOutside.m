function [c, ceq] = nTrabPointsOutside(x, shapeModel, outerShapePts, outerShapeF)

    %% This script serves as the function input to the fmincon procedure
    %  that provides the constraints for ensuring there are no points of
    %  the trabecular surface outside of the tibia
    
    %  Inputs:
    %       x - n x 1 array of PC scores to use in reconstruction coming from the optimisation
    %       shapeModel - the structure with shape model data for reconstructing the trabecular surface
    %       outerShapePts - n x 3 array of the XYZ outer tibia points to search for outer trabecular points
    %       outerShapeF - faces corresponding to the outer points for surface reconstruction
    
    %  Outputs:
    %       c - blank output as no inequality constraints
    %       ceq - equality constraints of the number of trabecular points outside (must equal zero)
    
    %% Set-up
    
    %Set blank inequality constraint
    c = [];
    
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
        
    %Identify trabecular points within (and hence also outside) of tibia
    trabPtsIn = intriangulation(outerShapePts, outerShapeF, finalPoints, 0);
    
    %Calculate number of points outside and set as equality constraint
    ceq = sum(~trabPtsIn);
    
end