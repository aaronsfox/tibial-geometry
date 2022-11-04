function [jaccardSimilarity] = calcJaccardShapeModel(actualV, predictedV, shapeModel, nParts, preSplit, noNodesPart1)

    %% Function to calculate the Jaccard index - i.e. a measure of volume
    %  consistency - between the actual surface and reconsturcted surface.
    %
    %  This achieved by creating points within a volume for each surfaces
    %  and contrasting to determining matches and calculate Jaccard index.    
        
    %  Inputs:
    %       actualV - n x 3 array of the XYZ points of the actual surface
    %       predictedV - n x 3 array of the XYZ points of the predicted surface
    %       shapeModel - the structure with shape model data for the surface 
    %       nParts - options of 1 or 2 for single bone vs. dual surfaces
    %       preSplitFaces - only relevant for 2 part option. True if faces are already split into F1 & F2 structures
    %       noNodesPart1 - only relevant for 2 part options that have split faces. Identifies number of nodes to split parts by.
    
    if nParts == 2
        
        if preSplit
            
            %Extract nodes and faces
            
            %Actual
            %Part 1
            actualF1 = shapeModel.F1;
            actualV1 = actualV(1:noNodesPart1,:);
            %Part 2
            actualF2 = shapeModel.F2;
            actualV2 = actualV(noNodesPart1+1:end,:);
            
            %Predicted
            %Part 1
            predictedF1 = shapeModel.F1;
            predictedV1 = predictedV(1:noNodesPart1,:);
            %Part 2
            predictedF2 = shapeModel.F2;
            predictedV2 = predictedV(noNodesPart1+1:end,:);            
            
        else
    
            %Group surfaces
            %Actual
            [groupIndexVertices,groupIndexFaces] = groupVertices(shapeModel.F, actualV,0);        
            [actualF1, actualV1] = patchCleanUnused(shapeModel.F(groupIndexFaces == 1,:), actualV);
            [actualF2, actualV2] = patchCleanUnused(shapeModel.F(groupIndexFaces == 2,:), actualV);
            %Predicted
            [groupIndexVertices,groupIndexFaces] = groupVertices(shapeModel.F, predictedV,0);        
            [predictedF1, predictedV1] = patchCleanUnused(shapeModel.F(groupIndexFaces == 1,:), predictedV);
            [predictedF2, predictedV2] = patchCleanUnused(shapeModel.F(groupIndexFaces == 2,:), predictedV);
            
        end

        %Set ggremesh options structure
        opts1.nb_pts = length(actualV1); %Set desired number of points
        opts1.disp_on = 0; % Turn off command window text display
        opts2.nb_pts = length(actualV2); %Set desired number of points
        opts2.disp_on = 0; % Turn off command window text display

        %Fix surfaces to ensure appropriate meshing
        %Actual
        [actualF1, actualV1] = mergeVertices(actualF1, actualV1);
        [actualF1, actualV1] = ggremesh(actualF1, actualV1, opts1);
        [actualF2, actualV2] = mergeVertices(actualF2, actualV2);
        [actualF2, actualV2] = ggremesh(actualF2, actualV2, opts2);
        %Predicted
        [predictedF1, predictedV1] = mergeVertices(predictedF1, predictedV1);
        [predictedF1, predictedV1] = ggremesh(predictedF1, predictedV1, opts1);
        [predictedF2, predictedV2] = mergeVertices(predictedF2, predictedV2);
        [predictedF2, predictedV2] = ggremesh(predictedF2, predictedV2, opts2);

        %Find interior point
        innerPointActual1 = getInnerPoint(actualF1, actualV1);
        innerPointActual2 = getInnerPoint(actualF2, actualV2);
        innerPointPredicted1 = getInnerPoint(predictedF1, predictedV1);
        innerPointPredicted2 = getInnerPoint(predictedF2, predictedV2);

        %Create volumetric meshes using tetGen

        %Original surface
        tetVolume = tetVolMeanEst(actualF1, actualV1); %Volume for regular tets
        tetGenStruct.stringOpt = '-pq1.2AaY';
        tetGenStruct.Faces = actualF1;
        tetGenStruct.Nodes = actualV1;
        tetGenStruct.holePoints = [];
        tetGenStruct.faceBoundaryMarker = ones(size(actualF1,1),1); %Face boundary markers
        tetGenStruct.regionPoints = innerPointActual1; %region points
        tetGenStruct.regionA = tetVolume;
        [meshActual1] = runTetGen(tetGenStruct); %Run tetGen

        tetVolume = tetVolMeanEst(actualF2, actualV2); %Volume for regular tets
        tetGenStruct.stringOpt = '-pq1.2AaY';
        tetGenStruct.Faces = actualF2;
        tetGenStruct.Nodes = actualV2;
        tetGenStruct.holePoints = [];
        tetGenStruct.faceBoundaryMarker = ones(size(actualF2,1),1); %Face boundary markers
        tetGenStruct.regionPoints = innerPointActual2; %region points
        tetGenStruct.regionA = tetVolume;
        [meshActual2] = runTetGen(tetGenStruct); %Run tetGen

        %Predicted surface
        tetVolume = tetVolMeanEst(predictedF1, predictedV1); %Volume for regular tets
        tetGenStruct.stringOpt = '-pq1.2AaY';
        tetGenStruct.Faces = predictedF1;
        tetGenStruct.Nodes = predictedV1;
        tetGenStruct.holePoints = [];
        tetGenStruct.faceBoundaryMarker = ones(size(predictedF1,1),1); %Face boundary markers
        tetGenStruct.regionPoints = innerPointPredicted1; %region points
        tetGenStruct.regionA = tetVolume;
        [meshPredicted1] = runTetGen(tetGenStruct); %Run tetGen

        tetVolume = tetVolMeanEst(predictedF2, predictedV2); %Volume for regular tets
        tetGenStruct.stringOpt = '-pq1.2AaY';
        tetGenStruct.Faces = predictedF2;
        tetGenStruct.Nodes = predictedV2;
        tetGenStruct.holePoints = [];
        tetGenStruct.faceBoundaryMarker = ones(size(predictedF2,1),1); %Face boundary markers
        tetGenStruct.regionPoints = innerPointPredicted2; %region points
        tetGenStruct.regionA = tetVolume;
        [meshPredicted2] = runTetGen(tetGenStruct); %Run tetGen

        %Identify number of reconstructed volume points that are inside original
        %i.e. true positives
        reconInOrig1 = intriangulation(meshActual1.nodes, meshActual1.facesBoundary, meshPredicted1.nodes, 0);
        TP1 = sum(reconInOrig1);

        %Identify number of reconstructed volume that are outside original
        %i.e. false positives
        FP1 = sum(~reconInOrig1);

        %Identify number of original points not in reconstructed volume
        %i.e. false negatives
        origInRecon2 = intriangulation(meshPredicted2.nodes, meshPredicted2.facesBoundary, meshActual2.nodes, 0);
        FN1 = sum(~origInRecon2);

        %i.e. true positives
        reconInOrig2 = intriangulation(meshActual2.nodes, meshActual2.facesBoundary, meshPredicted2.nodes, 0);
        TP2 = sum(reconInOrig2);

        %Identify number of reconstructed volume that are outside original
        %i.e. false positives
        FP2 = sum(~reconInOrig2);

        %Identify number of original points not in reconstructed volume
        %i.e. false negatives
        origInRecon2 = intriangulation(meshPredicted2.nodes, meshPredicted2.facesBoundary, meshActual2.nodes, 0);
        FN2 = sum(~origInRecon2);

        %Calculate jaccard similarity
        jaccardSimilarity = (TP1+TP2) / ((TP1+TP2) + (FP1+FP2) + (FN1+FN2));
        
    elseif nParts == 1    
    
        %% Single part option    

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
        
    else
        
        error('Only 1 or 2 part options are available.')
        
    end
    
end