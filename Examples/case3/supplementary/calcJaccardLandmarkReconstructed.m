function [jaccardSimilarityTibia, jaccardSimilarityFibula, jaccardSimilarityAll] = calcJaccardLandmarkReconstructed(tibiaF, fibulaF, tibiaV, fibulaV, tibiaReconstructedV, fibulaReconstructedV)

    %% Function to calculate the Jaccard index - i.e. a measure of volume
    %  consistency - between the actual surface and reconsturcted surface.
    %
    %  This achieved by creating points within a volume for each surfaces
    %  and contrasting to determining matches and calculate Jaccard index.    
        
    %  Inputs:
    %       tibiaF - surface faces relating to the tibia nodes
    %       fibulaF - surface faces relating to the fibula nodes
    %       tibiaV - n x 3 array of the XYZ points of the actual tibia surface 
    %       fibulaV - n x 3 array of the XYZ points of the actual fibula surface
    %       tibiaReconstructedV - n x 3 array of the XYZ points of the reconstructed tibia surface 
    %       fibulaReconstructedV - n x 3 array of the XYZ points of the reconstructed fibula surface
    
    %% Tibia    

    %Set ggremesh options structure
    opts.nb_pts = length(tibiaV); %Set desired number of points
    opts.disp_on = 0; % Turn off command window text display

    %Fix surfaces to ensure appropriate meshing
    %Actual
    [tibiaActualF, tibiaActualV] = mergeVertices(tibiaF, tibiaV);
    [tibiaActualF, tibiaActualV] = ggremesh(tibiaActualF, tibiaActualV, opts);
    %Predicted
    [tibiaReconstructedF, tibiaReconstructedV] = mergeVertices(tibiaF, tibiaReconstructedV);
    [tibiaReconstructedF, tibiaReconstructedV] = ggremesh(tibiaReconstructedF, tibiaReconstructedV, opts);

    %Find interior point
    innerPointActual = getInnerPoint(tibiaActualF, tibiaActualV);
    innerPointPredicted = getInnerPoint(tibiaReconstructedF, tibiaReconstructedV);

    %Create volumetric meshes using tetGen

    %Original surface
    tetVolume = tetVolMeanEst(tibiaActualF, tibiaActualV); %Volume for regular tets
    tetGenStruct.stringOpt = '-pq1.2AaY';
    tetGenStruct.Faces = tibiaActualF;
    tetGenStruct.Nodes = tibiaActualV;
    tetGenStruct.holePoints = [];
    tetGenStruct.faceBoundaryMarker = ones(size(tibiaActualF,1),1); %Face boundary markers
    tetGenStruct.regionPoints = innerPointActual; %region points
    tetGenStruct.regionA = tetVolume;
    [meshActual] = runTetGen(tetGenStruct); %Run tetGen

    %Reconstructed surface
    tetVolume = tetVolMeanEst(tibiaReconstructedF, tibiaReconstructedV); %Volume for regular tets
    tetGenStruct.stringOpt = '-pq1.2AaY';
    tetGenStruct.Faces = tibiaReconstructedF;
    tetGenStruct.Nodes = tibiaReconstructedV;
    tetGenStruct.holePoints = [];
    tetGenStruct.faceBoundaryMarker = ones(size(tibiaReconstructedF,1),1); %Face boundary markers
    tetGenStruct.regionPoints = innerPointPredicted; %region points
    tetGenStruct.regionA = tetVolume;
    [meshPredicted] = runTetGen(tetGenStruct); %Run tetGen

    %Identify number of reconstructed volume points that are inside original
    %i.e. true positives
    reconInOrig = intriangulation(meshActual.nodes, meshActual.facesBoundary, meshPredicted.nodes, 0);
    tibiaTP = sum(reconInOrig);

    %Identify number of reconstructed volume that are outside original
    %i.e. false positives
    tibiaFP = sum(~reconInOrig);

    %Identify number of original points not in reconstructed volume
    %i.e. false negatives
    origInRecon = intriangulation(meshPredicted.nodes, meshPredicted.facesBoundary, meshActual.nodes, 0);
    tibiaFN = sum(~origInRecon);

    %Calculate jaccard similarity
    jaccardSimilarityTibia = tibiaTP / (tibiaTP + tibiaFP + tibiaFN);
    
    %% Fibula
    
    %Reset ggremesh options structure
    opts.nb_pts = length(fibulaV); %Set desired number of points
    
    %Fix surfaces to ensure appropriate meshing
    %Actual
    [fibulaActualF, fibulaActualV] = mergeVertices(fibulaF, fibulaV);
    [fibulaActualF, fibulaActualV] = ggremesh(fibulaActualF, fibulaActualV, opts);
    %Predicted
    [fibulaReconstructedF, fibulaReconstructedV] = mergeVertices(fibulaF, fibulaReconstructedV);
    [fibulaReconstructedF, fibulaReconstructedV] = ggremesh(fibulaReconstructedF, fibulaReconstructedV, opts);

    %Find interior point
    innerPointActual = getInnerPoint(fibulaActualF, fibulaActualV);
    innerPointPredicted = getInnerPoint(fibulaReconstructedF, fibulaReconstructedV);

    %Create volumetric meshes using tetGen

    %Original surface
    tetVolume = tetVolMeanEst(fibulaActualF, fibulaActualV); %Volume for regular tets
    tetGenStruct.stringOpt = '-pq1.2AaY';
    tetGenStruct.Faces = fibulaActualF;
    tetGenStruct.Nodes = fibulaActualV;
    tetGenStruct.holePoints = [];
    tetGenStruct.faceBoundaryMarker = ones(size(fibulaActualF,1),1); %Face boundary markers
    tetGenStruct.regionPoints = innerPointActual; %region points
    tetGenStruct.regionA = tetVolume;
    [meshActual] = runTetGen(tetGenStruct); %Run tetGen

    %Reconstructed surface
    tetVolume = tetVolMeanEst(fibulaReconstructedF, fibulaReconstructedV); %Volume for regular tets
    tetGenStruct.stringOpt = '-pq1.2AaY';
    tetGenStruct.Faces = fibulaReconstructedF;
    tetGenStruct.Nodes = fibulaReconstructedV;
    tetGenStruct.holePoints = [];
    tetGenStruct.faceBoundaryMarker = ones(size(fibulaReconstructedF,1),1); %Face boundary markers
    tetGenStruct.regionPoints = innerPointPredicted; %region points
    tetGenStruct.regionA = tetVolume;
    [meshPredicted] = runTetGen(tetGenStruct); %Run tetGen

    %Identify number of reconstructed volume points that are inside original
    %i.e. true positives
    reconInOrig = intriangulation(meshActual.nodes, meshActual.facesBoundary, meshPredicted.nodes, 0);
    fibulaTP = sum(reconInOrig);

    %Identify number of reconstructed volume that are outside original
    %i.e. false positives
    fibulaFP = sum(~reconInOrig);

    %Identify number of original points not in reconstructed volume
    %i.e. false negatives
    origInRecon = intriangulation(meshPredicted.nodes, meshPredicted.facesBoundary, meshActual.nodes, 0);
    fibulaFN = sum(~origInRecon);

    %Calculate jaccard similarity
    jaccardSimilarityFibula = fibulaTP / (fibulaTP + fibulaFP + fibulaFN);
    
    %% Calculate overall Jaccard similarity
    
    jaccardSimilarityAll = (tibiaTP+fibulaTP) / ((tibiaTP+fibulaTP) + (tibiaFP+fibulaFP) + (tibiaFN+fibulaFN));
    
end