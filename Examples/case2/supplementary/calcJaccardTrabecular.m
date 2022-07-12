function [jaccardSimilarity] = calcJaccardTrabecular(actualV, predictedV, shapeModel)

    %% Function to calculate the Jaccard index - i.e. a measure of volume
    %  consistency - between the actual trabecular and predicted trabecular
    
        
    %  Inputs:
    %       actualV - n x 3 array of the XYZ points of the actual trabecular
    %       predictedV - n x 3 array of the XYZ points of the predicted trabecular
    %       shapeModel - the structure with shape model data for the trabecular    
    
    %% Calculate Jaccard similarity    

    %Set ggremesh options structure
    opts.nb_pts = length(shapeModel.meanPoints); %Set desired number of points
    opts.disp_on = 0; % Turn off command window text display

    %Fix surfaces to ensure appropriate meshing
    %Actual
    [actualF, actualV] = mergeVertices(shapeModel.F, actualV);
    [actualF, actualV] = ggremesh(actualF, actualV, opts);
    %Predicted
    [predictedF, predictedV] = mergeVertices(shapeModel.F, predictedV);
    [predictedF, predictedV] = ggremesh(predictedF, predictedV, opts);

    %Find interior point
    innerPointActual = getInnerPoint(actualF, actualV);
    innerPointPredicted = getInnerPoint(predictedF, predictedV);

    %Create volumetric meshes using tetGen

    %Original surface
    tetVolume = tetVolMeanEst(actualF, actualV); %Volume for regular tets
    tetGenStruct.stringOpt = '-pq1.2AaY';
    tetGenStruct.Faces = actualF;
    tetGenStruct.Nodes = actualV;
    tetGenStruct.holePoints = [];
    tetGenStruct.faceBoundaryMarker = ones(size(actualF,1),1); %Face boundary markers
    tetGenStruct.regionPoints = innerPointActual; %region points
    tetGenStruct.regionA = tetVolume;
    [meshActual] = runTetGen(tetGenStruct); %Run tetGen

    %Reconstructed surface
    tetVolume = tetVolMeanEst(predictedF, predictedV); %Volume for regular tets
    tetGenStruct.stringOpt = '-pq1.2AaY';
    tetGenStruct.Faces = predictedF;
    tetGenStruct.Nodes = predictedV;
    tetGenStruct.holePoints = [];
    tetGenStruct.faceBoundaryMarker = ones(size(predictedF,1),1); %Face boundary markers
    tetGenStruct.regionPoints = innerPointPredicted; %region points
    tetGenStruct.regionA = tetVolume;
    [meshPredicted] = runTetGen(tetGenStruct); %Run tetGen
    
    %Check if reconstructed surface meshes
    %Sometimes needs smoothing to sort out intersecting faces
    if isempty(meshPredicted.nodes)

        %Smoothing surface
        smoothPar.n = 2; %2 iterations
        smoothPar.Method = 'LAP'; %Laplacian smoothing method        
        smoothPredictedV = patchSmooth(predictedF, predictedV,[],smoothPar);

        %Remesh using ggremesh to reduce complexity
        [smoothPredictedF,smoothPredictedV] = ggremesh(predictedF, smoothPredictedV, opts);

        %Get an updated inner point estimate
        innerPointPredicted = getInnerPoint(smoothPredictedF,smoothPredictedV);

        %Mesh newly smoothed reconstructed surface
        tetVolume = tetVolMeanEst(smoothPredictedF,smoothPredictedV); %Volume for regular tets
        tetGenStruct.stringOpt = '-pq1.2AaY';
        tetGenStruct.Faces = smoothPredictedF;
        tetGenStruct.Nodes = smoothPredictedV;
        tetGenStruct.holePoints = [];
        tetGenStruct.faceBoundaryMarker = ones(size(smoothPredictedF,1),1); %Face boundary markers
        tetGenStruct.regionPoints = innerPointPredicted; %region points
        tetGenStruct.regionA = tetVolume;
        [meshPredicted] = runTetGen(tetGenStruct); %Run tetGen
    
    end

    %Identify number of reconstructed volume points that are inside original
    %i.e. true positives
    reconInOrig = intriangulation(meshActual.nodes, meshActual.facesBoundary, meshPredicted.nodes, 0);
    TP = sum(reconInOrig);

    %Identify number of reconstructed volume that are outside original
    %i.e. false positives
    FP = sum(~reconInOrig);

    %Identify number of original points not in reconstructed volume
    %i.e. false negatives
    origInRecon = intriangulation(meshPredicted.nodes, meshPredicted.facesBoundary, meshActual.nodes, 0);
    FN = sum(~origInRecon);

    %Calculate jaccard similarity
    jaccardSimilarity = TP / (TP + FP + FN);

end