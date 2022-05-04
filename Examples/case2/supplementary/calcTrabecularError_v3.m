function [sumErrors] = calcTrabecularError_v3(x, predictedTrabV, tibiaV, shapeModel, w)

    %% This script serves as the function input to the fmincon procedure
    %  where we aim to minimise the scale factors in correcting the
    %  trabecular surface while combining with an equality constraint of
    %  keeping all the points inside the tibia.
    %
    %  Effectively the aim is to minimise the scaling and translation
    %  applied to the surface to ensure that it fits within the tibia outer
    %  surface.
    
    %  Inputs:
    %       x - 6 x 1 array of:
    %               3 values related to scaling of surface
    %               3 values related to XYZ rotation of surface (in degrees)
    %               3 values related to XYZ translation of surface (in mm)
    %       predictedTrabV - n x 3 array of the XYZ points of the predictd trabecular
    %       tibiaV - n x 3 array of the XYZ points of the outer tibia surface 
    %       shapeModel - the structure with shape model data for the tibia
    %       w - weights applied to the error of the inputs
    %               3 weights for XYZ scaling parameter
    %               3 weights for XYZ translation parameter
    %               3 weights for XYZ rotations parameter
    %               1 weight for the number of points outside boundary
       
    %% Manipulate surface
    
    %Apply scaling
    trabPoints = x(1:3) .* predictedTrabV;

    %Recentre trabecular to tibia
    
    %Get central points along both point clouds
    centralTibia = mean([max(tibiaV);min(tibiaV)]);
    centralTrab = mean([max(trabPoints);min(trabPoints)]);    

    %Adjust the predicted trabecular to centre of tibia
    trabPoints = trabPoints - (centralTrab - centralTibia);
    
    %Get the error in scale factors relative to neutral (i.e. 1)
    error(1) = abs(1 - x(1)) * w(1);
    error(2) = abs(1 - x(2)) * w(2);
    error(3) = abs(1 - x(3)) * w(3);   
    
    %Apply rotations to the trabecular surface
    trabPoints = transformPoint3d(trabPoints, ...
        createRotationOx(deg2rad(x(4))));
    trabPoints = transformPoint3d(trabPoints, ...
        createRotationOy(deg2rad(x(5))));
    trabPoints = transformPoint3d(trabPoints, ...
        createRotationOz(deg2rad(x(6))));
    
    %Get the error in rotational factors relative to neutral (i.e. 0)
    error(4) = abs(0 - x(4)) * w(4);
    error(5) = abs(0 - x(5)) * w(5);
    error(6) = abs(0 - x(6)) * w(6);
    
    %Apply translation to the trabecular surface
    trabPoints = trabPoints - x(7:9);
    
    %Get the error in translational factors relative to neutral (i.e. 0)
    error(7) = abs(0 - x(7)) * w(7);
    error(8) = abs(0 - x(8)) * w(8);
    error(9) = abs(0 - x(9)) * w(9);
    
    %Determine number of points outside tibia boundary
    trabPtsIn = intriangulation(tibiaV, shapeModel.F, predictedTrabV, 0);
    error(10) = sum(~trabPtsIn) * w(10);
    
    %Sum the scale, translational and rotational errors + points outside
    sumErrors = sum(error);
    
end