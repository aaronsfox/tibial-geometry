function [c, ceq] = nTrabPointsOutside_v4(x, predictedTrabV, tibiaV, shapeModel)

    %% This script serves as the function input to the fmincon procedure
    %  that provides the constraints for ensuring there are no points of
    %  the trabecular surface outside of the tibia
    
    %  Inputs:
    %       x - 10500 data points corresponding to distance travelled in 3D by points 
    %       predictedTrabV - n x 3 array of the XYZ points of the predicted trabecular
    %       tibiaV - n x 3 array of the XYZ points of the outer tibia surface
    %       shapeModel - the structure with shape model data for the tibia
    
    %  Outputs:
    %       c - blank output as no inequality constraints
    %       ceq - equality constraints of the number of trabecular points outside (must equal zero)
    
    %% Set-up
    
    %Set blank inequality constraint
    c = [];
    
    %% Manipulate surface
    
    %Reshape the translation points to 3D
    transPts = reshape(x, [length(x)/3,3]);
    
    %Apply translation
    newV = predictedTrabV - transPts;
    
    %% Check points inside tibia triangulation
        
    %Identify trabecular points within (and hence also outside) of tibia
    trabPtsIn = intriangulation(tibiaV, shapeModel.F, newV, 0);
    
    %Calculate number of points outside and set as equality constraint
    ceq = sum(~trabPtsIn);
    
end