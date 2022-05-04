function [c, ceq] = nTrabPointsOutside_v2(x, origPoints, outerShapePts, outerShapeF)

    %% This script serves as the function input to the fmincon procedure
    %  that provides the constraints for ensuring there are no points of
    %  the trabecular surface outside of the tibia
    
    %  Inputs:
    %       x - n x 1 array of values to use in reconstruction coming from the optimisation
    %       shapeModel - the structure with shape model data for reconstructing the trabecular surface
    %       outerShapePts - n x 3 array of the XYZ outer tibia points to search for outer trabecular points
    %       outerShapeF - faces corresponding to the outer points for surface reconstruction
    
    %  Outputs:
    %       c - blank output as no inequality constraints
    %       ceq - equality constraints of the number of trabecular points outside (must equal zero)
    
    %% Set-up
    
    %Set blank inequality constraint
    c = [];

    %Extract the scaling factor
    scaleFactor(1:3,1) = x(1:3);
    
    %Extract the translational matrix
    transMat(1:3,1) = x(4:6);
    
    %Extract the rotation matrix
    rotMat(1,1:3) = x(7:9);
    rotMat(2,1:3) = x(10:12);
    rotMat(3,1:3) = x(13:15);
    
    %% Reconstruct surface
    
    %Apply translation, rotation and scaling
    finalPoints = scaleFactor' .* origPoints * rotMat' + ...
        repmat(transMat',[length(origPoints) 1]);
        
    %Identify trabecular points within (and hence also outside) of tibia
    trabPtsIn = intriangulation(outerShapePts, outerShapeF, finalPoints, 0);
    
    %Calculate number of points outside and set as equality constraint
    ceq = sum(~trabPtsIn);
    
end