function [sumPointDistance] = calcTrabecularError_v4(x, predictedTrabV)

    %% This script serves as the function input to the fmincon procedure
    %  where we aim to minimise the distance travelled by points of the
    %  trabecular while constraining them to be within the boundaries of
    %  the tibia.
    
    %  Inputs:
    %       x - 10500 data points corresponding to distance travelled in 3D by points 
    %       predictedTrabV - n x 3 array of the XYZ points of the predicted trabecular
       
    %% Manipulate surface
    
    %Reshape the translation points to 3D
    transPts = reshape(x, [length(x)/3,3]);
    
    %Apply translation
    newV = predictedTrabV - transPts;
    
    %Calculate the distance travelled by each point
    pointDistance = distancePoints3d(predictedTrabV, newV);
    
    %Sum the point distances as the value to minimise and convert to m
    sumPointDistance = sum(pointDistance) / 1000;
    
end