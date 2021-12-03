function calcRelativeTrabecularVol(tibiaF, tibiaV, trabF, trabV)

    %% Function to clean up main code where we calculate the relative trabecular
    %  volume of the tibia along its length at 5% intervals
    
    %Determine max and min points of tibia surface along Y-axis
    maxTibLength = max(tibiaV(:,2));
    minTibLength = min(tibiaV(:,2));

    %Determine 0-100% at 5% intervals along Y-axis tibia length
    tibLengthPts = linspace(minTibLength, maxTibLength, 100/5+1);

    %Loop through the pairs of points along tibia length
    for pt = 1:length(tibLengthPts)-1

        %Get the two current points
        minLevel = tibLengthPts(pt);
        maxLevel = tibLengthPts(pt+1);

        %Create the logic for cutting away faces

        %Distal cut

        %Cortical
        logicVerticesCortical = tibiaV(:,2) >= minLevel; %%% & (tibiaV(:,2) <= maxLevel);
        logicFacesCortical = all(logicVerticesCortical(tibiaF),2);
        logicFacesCortical = triSurfLogicSharpFix(tibiaF, logicFacesCortical, 3); %alter logic for smoothness
        %Trabecular
        logicVerticesTrabecular = (trabV(:,2) >= minLevel); %%% & (trabV(:,2) <= maxLevel);
        logicFacesTrabecular = all(logicVerticesTrabecular(trabF),2);
        logicFacesTrabecular = triSurfLogicSharpFix(trabF, logicFacesTrabecular, 3); %alter logic for smoothness

% % %             %Visualise
% % %             cFigure; hold on;
% % %             gpatch(tibiaF, tibiaV, logicFacesCortical, 'none', 0.5);
% % %             gpatch(trabF, trabV, logicFacesTrabecular, 'k');
% % %             camlight('headlight');
% % %             axisGeom();

        %Cut away the faces using the logic determined above
        %Faces to keep
        tibiaF = tibiaF(logicFacesCortical,:);
        trabF = trabF(logicFacesTrabecular,:);
        %Remove unused points
        [tibiaF, tibiaV] = patchCleanUnused(tibiaF, tibiaV);
        [trabF, trabV] = patchCleanUnused(trabF, trabV);

        %Attempt to self triangulate potentially jagged edge
        %Cortical
        tibiaEb = patchBoundary(tibiaF, tibiaV); %Get boundary edges
        tibiaIndBoundary = edgeListToCurve(tibiaEb); %Convert boundary edges to a curve list
        tibiaIndBoundary = tibiaIndBoundary(1:end-1); %Trim off last point since it is equal to first on a closed loop
        angleThreshold = pi*(120/180); %threshold for self triangulation
        [tibiaF, tibiaV, indBoundaryCortical] = ...
            triSurfSelfTriangulateBoundary(tibiaF, tibiaV, tibiaIndBoundary, angleThreshold, 1);
        %Trabecular
        trabEb = patchBoundary(trabF, trabV); %Get boundary edges
        trabIndBoundary = edgeListToCurve(trabEb); %Convert boundary edges to a curve list
        trabIndBoundary = trabIndBoundary(1:end-1); %Trim off last point since it is equal to first on a closed loop
        angleThreshold = pi*(120/180); %threshold for self triangulation
        [trabF, trabV, indBoundaryTrabecular] = ...
            triSurfSelfTriangulateBoundary(trabF, trabV, trabIndBoundary, angleThreshold, 1);

        %Convert edge points to have the same y-edge level as the cut
        tibiaV(indBoundaryCortical,2) = minLevel;
        trabV(indBoundaryTrabecular,2) = minLevel;

% % %             %Visualise
% % %             cFigure; hold on;
% % %             %Cortical
% % %             gpatch(keepTibiaF, keepTibiaV, 'gw','none', 0.3);
% % %             plotV(keepTibiaV(indBoundaryCortical,:),'r-','LineWidth', 1.5);
% % %             %Trabecular
% % %             gpatch(keepTrabF, keepTrabV, 'rw','k', 1);
% % %             plotV(keepTrabV(indBoundaryTrabecular,:),'b-','LineWidth', 1.5);
% % %             axisGeom();            
% % %             camlight('headlight');

        %Close over distal regions
        %Cortical
        %Remesh in 2D with X and Z axes
        [distalF, distalV] = regionTriMesh2D({tibiaV(indBoundaryCortical, [1,3])}, ...
            mean(patchEdgeLengths(tibiaF, tibiaV)) / 2, 0, 0);
        %Fix up the axes so that it aligns with original surface
        distalV(:,1:3) = [distalV(:,1), ones(length(distalV),1)*minLevel, distalV(:,2)];
        %Trabecular
        %Remesh in 2D with X and Z axes
        [distalTrabF, distalTrabV] = regionTriMesh2D({trabV(indBoundaryTrabecular, [1,3])}, ...
            mean(patchEdgeLengths(trabF, trabV)) / 2, 0, 0);
        %Fix up the axes so that it aligns with original surface
        distalTrabV(:,1:3) = [distalTrabV(:,1), ones(length(distalTrabV),1)*minLevel, distalTrabV(:,2)];

% % %         %Visualise
% % %         cFigure; hold on;
% % %         %Cortical
% % %         gpatch(tibiaF, tibiaV, 'gw','none', 0.3);
% % %         plotV(tibiaV(indBoundaryCortical,:),'r-','LineWidth', 1.5);
% % %         gpatch(distalF, distalV, 'gw','none', 0.3);
% % %         %Trabecular
% % %         gpatch(trabF, trabV, 'rw','k', 1);
% % %         plotV(trabV(indBoundaryTrabecular,:),'b-','LineWidth', 1.5);
% % %         gpatch(distalTrabF, distalTrabV, 'rw','k', 1);
% % %         axisGeom();            
% % %         camlight('headlight');

        %Merge distal ends with main surfaces
        [tibiaF, tibiaV] = joinElementSets({tibiaF, distalF},{tibiaV, distalV});
        [trabF, trabV] = joinElementSets({trabF, distalTrabF},{trabV, distalTrabV});
        
% % %             %Visualise
% % %             cFigure; hold on;
% % %             gpatch(tibiaF, tibiaV, 'gw', 'none', 0.3);
% % %             gpatch(trabF, trabV, 'rw', 'k');
% % %             camlight('headlight');
% % %             axisGeom();

        %Proximal cut
        
        %Cortical
        logicVerticesCortical2 = tibiaV(:,2) <= maxLevel;
        logicFacesCortical2 = all(logicVerticesCortical2(tibiaF),2);
        logicFacesCortical2 = triSurfLogicSharpFix(tibiaF, logicFacesCortical2, 3); %alter logic for smoothness
        %Trabecular
        logicVerticesTrabecular2 = trabV(:,2) <= maxLevel;
        logicFacesTrabecular2 = all(logicVerticesTrabecular2(trabF),2);
        logicFacesTrabecular2 = triSurfLogicSharpFix(trabF, logicFacesTrabecular2, 3); %alter logic for smoothness
        
        
        %%%%% TODO: SHELVING THIS FOR HERE - AS I DON'T THINK VOLUME EVERY
        %%%%% 5% IS APPROPIRATE, IT'S MORE LIKE THE THICKNESS AT EACH SORT
        %%%%% OF STEP UP THE TIBIA IS APPROPRIATE - I.E. CUT THE TIBIA
        %%%%% EVERY 1% AND ASSESS THE CORTICAL THICKNESS AT THAT CUT
        %%%%% LEVEL...?
        
        
% % %             %Visualise
% % %             cFigure; hold on;
% % %             gpatch(tibiaF, tibiaV, logicFacesCortical2, 'none', 0.5);
% % %             gpatch(trabF, trabV, logicFacesTrabecular2, 'k');
% % %             camlight('headlight');
% % %             axisGeom();

% % %         %Cut away the faces using the logic determined above
% % %         %Faces to keep
% % %         tibiaF = tibiaF(logicFacesCortical2,:);
% % %         trabF = trabF(logicFacesTrabecular2,:);
% % %         %Remove unused points
% % %         [tibiaF, tibiaV] = patchCleanUnused(tibiaF, tibiaV);
% % %         [trabF, trabV] = patchCleanUnused(trabF, trabV);
        
% % %             %Visualise
% % %             cFigure; hold on;
% % %             gpatch(tibiaF, tibiaV, 'gw', 'none', 0.3);
% % %             gpatch(trabF, trabV, 'rw', 'k');
% % %             camlight('headlight');
% % %             axisGeom();

% % %         %Attempt to self triangulate potentially jagged edge
% % %         %Cortical
% % %         tibiaEb2 = patchBoundary(tibiaF, tibiaV); %Get boundary edges
% % %         tibiaIndBoundary2 = edgeListToCurve(tibiaEb2); %Convert boundary edges to a curve list
% % %         tibiaIndBoundary = tibiaIndBoundary(1:end-1); %Trim off last point since it is equal to first on a closed loop
% % %         angleThreshold = pi*(120/180); %threshold for self triangulation
% % %         [tibiaF, tibiaV, indBoundaryCortical] = ...
% % %             triSurfSelfTriangulateBoundary(tibiaF, tibiaV, tibiaIndBoundary, angleThreshold, 1);
% % %         %Trabecular
% % %         trabEb = patchBoundary(trabF, trabV); %Get boundary edges
% % %         trabIndBoundary = edgeListToCurve(trabEb); %Convert boundary edges to a curve list
% % %         trabIndBoundary = trabIndBoundary(1:end-1); %Trim off last point since it is equal to first on a closed loop
% % %         angleThreshold = pi*(120/180); %threshold for self triangulation
% % %         [trabF, trabV, indBoundaryTrabecular] = ...
% % %             triSurfSelfTriangulateBoundary(trabF, trabV, trabIndBoundary, angleThreshold, 1);
% % % 
% % %         %Convert edge points to have the same y-edge level as the cut
% % %         tibiaV(indBoundaryCortical,2) = minLevel;
% % %         trabV(indBoundaryTrabecular,2) = minLevel;
% % % 
% % % % % %             %Visualise
% % % % % %             cFigure; hold on;
% % % % % %             %Cortical
% % % % % %             gpatch(keepTibiaF, keepTibiaV, 'gw','none', 0.3);
% % % % % %             plotV(keepTibiaV(indBoundaryCortical,:),'r-','LineWidth', 1.5);
% % % % % %             %Trabecular
% % % % % %             gpatch(keepTrabF, keepTrabV, 'rw','k', 1);
% % % % % %             plotV(keepTrabV(indBoundaryTrabecular,:),'b-','LineWidth', 1.5);
% % % % % %             axisGeom();            
% % % % % %             camlight('headlight');
% % % 
% % %         %Close over distal regions
% % %         %Cortical
% % %         %Remesh in 2D with X and Z axes
% % %         [distalF, distalV] = regionTriMesh2D({tibiaV(indBoundaryCortical, [1,3])}, ...
% % %             mean(patchEdgeLengths(tibiaF, tibiaV)) / 2, 0, 0);
% % %         %Fix up the axes so that it aligns with original surface
% % %         distalV(:,1:3) = [distalV(:,1), ones(length(distalV),1)*minLevel, distalV(:,2)];
% % %         %Trabecular
% % %         %Remesh in 2D with X and Z axes
% % %         [distalTrabF, distalTrabV] = regionTriMesh2D({trabV(indBoundaryTrabecular, [1,3])}, ...
% % %             mean(patchEdgeLengths(trabF, trabV)) / 2, 0, 0);
% % %         %Fix up the axes so that it aligns with original surface
% % %         distalTrabV(:,1:3) = [distalTrabV(:,1), ones(length(distalTrabV),1)*minLevel, distalTrabV(:,2)];

% % %         %Visualise
% % %         cFigure; hold on;
% % %         %Cortical
% % %         gpatch(tibiaF, tibiaV, 'gw','none', 0.3);
% % %         plotV(tibiaV(indBoundaryCortical,:),'r-','LineWidth', 1.5);
% % %         gpatch(distalF, distalV, 'gw','none', 0.3);
% % %         %Trabecular
% % %         gpatch(trabF, trabV, 'rw','k', 1);
% % %         plotV(trabV(indBoundaryTrabecular,:),'b-','LineWidth', 1.5);
% % %         gpatch(distalTrabF, distalTrabV, 'rw','k', 1);
% % %         axisGeom();            
% % %         camlight('headlight');
        
        
    end


end