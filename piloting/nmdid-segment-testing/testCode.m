
%Set path to a data folder
dataDir = '..\..\Segmentation\case-102480';
cd(dataDir);

%Load in a surface file
%Cortical tibia
[tibiaSTLstruct] = import_STL('102480-tibia-cortical.stl');
corticalF = tibiaSTLstruct.solidFaces{1}; %Faces
corticalV = tibiaSTLstruct.solidVertices{1}; %Vertices
[corticalF,corticalV] = mergeVertices(corticalF,corticalV);
%Trabecular tibia
[tibiaSTLstruct] = import_STL('102480-tibia-trabecular.stl');
trabecularF = tibiaSTLstruct.solidFaces{1}; %Faces
trabecularV = tibiaSTLstruct.solidVertices{1}; %Vertices
[trabecularF,trabecularV] = mergeVertices(trabecularF,trabecularV);
%Tibial plateau
[tibiaSTLstruct] = import_STL('102480-tibial-plateau.stl');
plateauF = tibiaSTLstruct.solidFaces{1}; %Faces
plateauV = tibiaSTLstruct.solidVertices{1}; %Vertices
[plateauF,plateauV] = mergeVertices(plateauF,plateauV);
%Ankle surface
[tibiaSTLstruct] = import_STL('102480-ankle-joint-surface.stl');
surfaceF = tibiaSTLstruct.solidFaces{1}; %Faces
surfaceV = tibiaSTLstruct.solidVertices{1}; %Vertices
[surfaceF,surfaceV] = mergeVertices(surfaceF,surfaceV);

%Check for holes in surfaces
[Eb] = patchBoundary(corticalF,corticalV);
if ~isempty(Eb)
    error('Holes in cortical surface')
end
[Eb] = patchBoundary(trabecularF,trabecularV);
if ~isempty(Eb)
    error('Holes in trabecular surface')
end

%Import malleoli
%xml_read looks like a custom function
tree = xml_read('malleoli.xml');
%%%% TODO: figure out ordering of landmarks otherwise it always needs to be
%%%% the same...
landmarks.MM = tree.Point(1).Coordinate;
landmarks.LM = tree.Point(2).Coordinate;
tree = xml_read('condyles.xml');
landmarks.LC = tree.Point(1).Coordinate;
landmarks.MC = tree.Point(2).Coordinate;
% % % tibPoints = [{'LC'},{'LM'},{'MC'},{'MM'}];
% % % for pp = 1:length(tibPoints)
% % %     tree = xml_read([tibPoints{pp},'.txt']);
% % %     landmarks.(char(tree.Point.Name)) = tree.Point.Coordinate;% / 1000;
% % %     clear tree
% % % end
% % % clear pp

%Create new landmarks
landmarks.IM = midPoint3d(landmarks.LM,landmarks.MM);
landmarks.IC = midPoint3d(landmarks.LC,landmarks.MC);

%% Visualise

cFigure; hold on
gpatch(trabecularF,trabecularV,'rw','none',0.5);
gpatch(corticalF,corticalV,'gw','none',0.5);
plotV(landmarks.LM,'r.');
plotV(landmarks.MM,'r.');
plotV(landmarks.LC,'r.');
plotV(landmarks.MC,'r.');
plotV(landmarks.IC,'r.');
plotV(landmarks.IM,'r.');
axisGeom; camlight headlight
title('Imported Surfaces and Landmarks');

%% Create and align tibial coordinate system

%Create the tibial planes
%Frontal
planes.frontal = createPlane(landmarks.IM,...
    landmarks.LC,landmarks.MC);
%Torsional
planes.torsional = createPlane(landmarks.IC,...
    landmarks.MM,landmarks.LM);

% % % drawPlane3d(planes.torsional);
% % % drawPlane3d(planes.frontal,'g');

%Create transform to get tibia aligned to the global plane
globalTransform = createBasisTransform3d('global',planes.torsional);

%Transform surfaces
for pp = 1:length(corticalV)
    corticalV(pp,:) = transformPoint3d(corticalV(pp,:),globalTransform);    
end
clear pp
for pp = 1:length(trabecularV)
    trabecularV(pp,:) = transformPoint3d(trabecularV(pp,:),globalTransform);    
end
clear pp
for pp = 1:length(plateauV)
    plateauV(pp,:) = transformPoint3d(plateauV(pp,:),globalTransform);    
end
clear pp
for pp = 1:length(surfaceV)
    surfaceV(pp,:) = transformPoint3d(surfaceV(pp,:),globalTransform);    
end
clear pp

%Transform landmarks
currLandmarks = fieldnames(landmarks);
for ff = 1:length(currLandmarks)
    landmarks.(currLandmarks{ff}) = transformPoint3d(landmarks.(currLandmarks{ff}),globalTransform);
end
clear ff

%Do the secondary rotation around the X-axis to make the tibia vertical
%along the Y-axis

%Identify distance between IC and IM, along XY plane
pt1 = [landmarks.IC(1),landmarks.IC(2)];
pt2 = [landmarks.IM(1),landmarks.IM(2)];
dIC_IM = sqrt((pt2(2) - pt1(2))^2 + (pt2(1) - pt1(1))^2);

%Identify distance of IC from IM along the x-axis
pt3 = [landmarks.IC(1),landmarks.IM(2)];
dIC_IMx = sqrt((pt3(2) - pt1(2))^2 + (pt3(1) - pt1(1))^2);

%Calculate angle to rotate about x-axis
rotAng = asin(dIC_IMx/dIC_IM);

%Create rotation matrix around Z-axis by specified angle (in radians)
rotZ = createRotationOz(rotAng*-1); %-ve for anti-clockwise

%Transform surfaces
for pp = 1:length(corticalV)
    corticalV(pp,:) = transformPoint3d(corticalV(pp,:),rotZ);    
end
clear pp
for pp = 1:length(trabecularV)
    trabecularV(pp,:) = transformPoint3d(trabecularV(pp,:),rotZ);    
end
clear pp
for pp = 1:length(plateauV)
    plateauV(pp,:) = transformPoint3d(plateauV(pp,:),rotZ);    
end
clear pp
for pp = 1:length(surfaceV)
    surfaceV(pp,:) = transformPoint3d(surfaceV(pp,:),rotZ);    
end
clear pp

%Transform landmarks
currLandmarks = fieldnames(landmarks);
for ff = 1:length(currLandmarks)
    landmarks.(currLandmarks{ff}) = transformPoint3d(landmarks.(currLandmarks{ff}),rotZ);
end
clear ff

%Create the transform to make the IM landmark the global origin
transMatrix = createTranslation3d([0,0,0] - landmarks.IM);

%Transform surfaces
for pp = 1:length(corticalV)
    corticalV(pp,:) = transformPoint3d(corticalV(pp,:),transMatrix);    
end
clear pp
for pp = 1:length(trabecularV)
    trabecularV(pp,:) = transformPoint3d(trabecularV(pp,:),transMatrix);    
end
clear pp
for pp = 1:length(plateauV)
    plateauV(pp,:) = transformPoint3d(plateauV(pp,:),transMatrix);    
end
clear pp
for pp = 1:length(surfaceV)
    surfaceV(pp,:) = transformPoint3d(surfaceV(pp,:),transMatrix);    
end
clear pp

%Transform landmarks
currLandmarks = fieldnames(landmarks);
for ff = 1:length(currLandmarks)
    landmarks.(currLandmarks{ff}) = transformPoint3d(landmarks.(currLandmarks{ff}),transMatrix);
end
clear ff

%The current bodies are aligned so that X is vertical and Y is lateral
%This needs to be shifted to align with the ISB recommendations so that Y
%is vertical and Z is lateral. This can easily be done by a few rotations
%about specific axes

%First, rotate about the z-axis by -90 degrees

%Transform surfaces
for pp = 1:length(corticalV)
    corticalV(pp,:) = transformPoint3d(corticalV(pp,:),createRotationOz(deg2rad(-90)));    
end
clear pp
for pp = 1:length(trabecularV)
    trabecularV(pp,:) = transformPoint3d(trabecularV(pp,:),createRotationOz(deg2rad(-90)));    
end
clear pp
for pp = 1:length(plateauV)
    plateauV(pp,:) = transformPoint3d(plateauV(pp,:),createRotationOz(deg2rad(-90)));    
end
clear pp
for pp = 1:length(surfaceV)
    surfaceV(pp,:) = transformPoint3d(surfaceV(pp,:),createRotationOz(deg2rad(-90)));    
end
clear pp

%Transform landmarks
currLandmarks = fieldnames(landmarks);
for ff = 1:length(currLandmarks)
    landmarks.(currLandmarks{ff}) = transformPoint3d(landmarks.(currLandmarks{ff}),createRotationOz(deg2rad(-90)));
end
clear ff

%Second, rotate about the y-axis by -90 degrees

%Transform surfaces
for pp = 1:length(corticalV)
    corticalV(pp,:) = transformPoint3d(corticalV(pp,:),createRotationOy(deg2rad(-90)));    
end
clear pp
for pp = 1:length(trabecularV)
    trabecularV(pp,:) = transformPoint3d(trabecularV(pp,:),createRotationOy(deg2rad(-90)));    
end
clear pp
for pp = 1:length(plateauV)
    plateauV(pp,:) = transformPoint3d(plateauV(pp,:),createRotationOy(deg2rad(-90)));    
end
clear pp
for pp = 1:length(surfaceV)
    surfaceV(pp,:) = transformPoint3d(surfaceV(pp,:),createRotationOy(deg2rad(-90)));    
end
clear pp

%Transform landmarks
currLandmarks = fieldnames(landmarks);
for ff = 1:length(currLandmarks)
    landmarks.(currLandmarks{ff}) = transformPoint3d(landmarks.(currLandmarks{ff}),createRotationOy(deg2rad(-90)));
end
clear ff

%% Revisualise

cFigure; hold on
gpatch(trabecularF,trabecularV,'rw','none',0.5);
gpatch(corticalF,corticalV,'gw','none',0.5);
plotV(landmarks.LM,'r.');
plotV(landmarks.MM,'r.');
plotV(landmarks.LC,'r.');
plotV(landmarks.MC,'r.');
plotV(landmarks.IC,'r.');
plotV(landmarks.IM,'r.');
axisGeom; camlight headlight
title('Imported Surfaces and Landmarks');

%% Get volumetric mesh of tibia

%From Example 3 in runTetGen

%%%%% IMPORTANT: for multi-region meshing, the order in joining the element
%%%%% sets and regions must be from outside-inside (i.e.
%%%%% cortical-trabecular) --- otherwise it won't work...

%Join surfaces
[F,V,C] = joinElementSets({corticalF trabecularF},{corticalV trabecularV});

%Mesh surfaces
[V_region1] = getInnerPoint(F(C==2,:),V); %First interior point
[V_region2] = getInnerPoint({corticalF trabecularF},{corticalV trabecularV}); %Second interior point
V_regions = [V_region1; V_region2]; %Collect region points
V_holes = []; %Define hole points
[regionTetVolume1] = tetVolMeanEst(corticalF,corticalV); %Volume estimate for regular tets
[regionTetVolume2] = tetVolMeanEst(trabecularF,trabecularV); %Volume estimate for regular tets
regionTetVolumes = [regionTetVolume1 regionTetVolume2];
stringOpt = '-pq1.2AaY'; %Tetgen options

%Create tetgen input structure
inputStruct.stringOpt = stringOpt; %Tetgen options
inputStruct.Faces = F; %Boundary faces
inputStruct.Nodes = V; %Nodes of boundary
inputStruct.faceBoundaryMarker = C;
inputStruct.regionPoints = V_regions; %Interior points for regions
inputStruct.holePoints = V_holes; %Interior points for holes
inputStruct.regionA = regionTetVolumes; %Desired tetrahedral volume for each region

% Mesh model using tetrahedral elements using tetGen
[meshOutput] = runTetGen(inputStruct); %Run tetGen

%Access mesh output structure
E = meshOutput.elements; %The elements
V = meshOutput.nodes; %The vertices or nodes
CE = meshOutput.elementMaterialID; %Element material or region id
Fb = meshOutput.facesBoundary; %The boundary faces
Cb = meshOutput.boundaryMarker; %The boundary markers

%Visualise
cFigure; hold on;
title('Input boundaries','FontSize',15);
hp(1) = gpatch(Fb,V,Cb,'none',0.3);
hp(2) = plotV(V_regions,'r.','MarkerSize',25);
legend(hp,{'Input mesh','Interior point(s)'},'Location','NorthWestOutside');
axisGeom(gca,15); camlight headlight;
colormap(cMap); icolorbar;

% Visualizing using |meshView|
meshView(meshOutput);







