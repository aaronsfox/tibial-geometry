function [c, ceq] = nTrabPointsOutside_v3(x, predictedTrabV, tibiaV, shapeModel)

    %% This script serves as the function input to the fmincon procedure
    %  that provides the constraints for ensuring there are no points of
    %  the trabecular surface outside of the tibia
    
    %  Inputs:
    %       x - 6 x 1 array of:
    %               3 values related to scaling of surface
    %               3 values related to XYZ rotation of surface (in degrees)
    %               3 values related to XYZ translation of surface (in mm)
    %       predictedTrabV - n x 3 array of the XYZ points of the predictd trabecular
    %       tibiaV - n x 3 array of the XYZ points of the outer tibia surface
    %       shapeModel - the structure with shape model data for the tibia
    
    %  Outputs:
    %       c - blank output as no inequality constraints
    %       ceq - equality constraints of the number of trabecular points outside (must equal zero)
    
    %% Set-up
    
    %Set blank inequality constraint
    c = [];
    
    %% Manipulate surface
    
    %Apply translation, rotation and scaling
    trabPoints = x(1:3) .* predictedTrabV;

    %Recentre trabecular to tibia
    
    %Get central points along both point clouds
    centralTibia = mean([max(tibiaV);min(tibiaV)]);
    centralTrab = mean([max(trabPoints);min(trabPoints)]);  
 
    %Adjust the predicted trabecular to centre of tibia
    trabPoints = trabPoints - (centralTrab - centralTibia);
    
    %Apply rotations to the trabecular surface
    trabPoints = transformPoint3d(trabPoints, ...
        createRotationOx(deg2rad(x(4))));
    trabPoints = transformPoint3d(trabPoints, ...
        createRotationOy(deg2rad(x(5))));
    trabPoints = transformPoint3d(trabPoints, ...
        createRotationOz(deg2rad(x(6))));
    
    %Apply translation to the trabecular surface
    trabPoints = trabPoints - x(7:9);
    
    %% Check points inside tibia triangulation
        
    %Identify trabecular points within (and hence also outside) of tibia
    trabPtsIn = intriangulation(tibiaV, shapeModel.F, trabPoints, 0);
    
    %Calculate number of points outside and set as equality constraint
    ceq = sum(~trabPtsIn);
    
end