%% This script provides a demonstration of how to generate the tibia-fibula
%  surface from a simple set of landmarks and the statistical shape model,
%  and provides an evaluation of the error associated with this across a
%  set of cases.
%
%  See the README.MD in this folder for more descriptive details on this
%  process.
%
% Authors:
%     Meghan Keast
%     Centre for Sport Research
%     Deakin University
%     mfkeast@deakin.edu.au
%
%     Aaron Fox
%     Centre for Sport Research
%     Deakin University
%     aaron.f@deakin.edu.au

%% Set-up

%Set home directory
homeDir = pwd;

%Add supplementary code paths
addpath(genpath('supplementary'));

%Load the tibia and trabecular shape models required
load('..\..\ShapeModels\tibia-fibula\tibiaFibulaShapeModel.mat');

% % % %Visualise mean shape model
% % % cFigure; hold on
% % % gpatch(tibiaFibulaShapeModel.F,tibiaFibulaShapeModel.meanPoints,'gw','k')
% % % axisGeom; camlight headlight

%Set options for remeshing
optionStruct_tib.nb_pts = 3500; %Set desired number of points
optionStruct_tib.disp_on = 0; % Turn off command window text display
optionStruct_fib.nb_pts = 2000; %Set desired number of points
optionStruct_fib.disp_on = 0; % Turn off command window text display

%% Load and prepare experimental data

%Set participant ID
pID = 'C01RTF';

%Tibia-fibula
[tibFibSTLstruct] = import_STL(['Nolte2016\',pID,'.stl']);
tibFibF = tibFibSTLstruct.solidFaces{1}; %Faces
tibFibV = tibFibSTLstruct.solidVertices{1}; %Vertices
[tibFibF,tibFibV] = mergeVertices(tibFibF,tibFibV);

%Split the surfaces to grab the tibia only

%Group the vertices and faces from the combined shape model
[groupIndexVertices,groupIndexFaces] = groupVertices(tibFibF, tibFibV,0);

%Identify which grouped section contains a higher volume
%This will be indicative of the tibia
if tetVolMeanEst(tibFibF(groupIndexFaces == 1,:),tibFibV) > ...
        tetVolMeanEst(tibFibF(groupIndexFaces == 2,:),tibFibV)
    %First index is tibia
    logicKeep = groupIndexFaces == 1;
else
    %Second index is tibia
    logicKeep = groupIndexFaces == 2;
end
%Separate the surfaces
[tibiaF, tibiaV] = patchCleanUnused(tibFibF(logicKeep,:), tibFibV);
[fibulaF, fibulaV] = patchCleanUnused(tibFibF(~logicKeep,:), tibFibV);

% % % %Visualise
% % % cFigure; hold on
% % % gpatch(tibiaF,tibiaV,'gw', 'none')
% % % gpatch(fibulaF,fibulaV,'bw', 'none')
% % % axisGeom; camlight headlight

%Extract the desired markers - RLKN, RMKN, RLM, RMM
%%%%% TODO: appropriate way to do this with Nolte dataset
sampleLandmarks.MC = [-100.4, -44.92, -446.5];
sampleLandmarks.LC = [-169.0, -18.91, -442.8];
sampleLandmarks.FH = [-158.9, 12.2, -467.2];
sampleLandmarks.MM = [-109, -52.16, -826.9];
sampleLandmarks.LM = [-140.2, 2.46, -840.9];
sampleLandmarks.TT = [-160.8, -61.22, -480.9];

%Create new landmarks needed to align tibia to global
sampleLandmarks.IM = midPoint3d(sampleLandmarks.LM,sampleLandmarks.MM);
sampleLandmarks.IC = midPoint3d(sampleLandmarks.LC,sampleLandmarks.MC);

%Remesh surfaces

%Run initial smoothing on surfaces
smoothPar.n = 2; %2 iterations
smoothPar.Method = 'LAP'; %Laplacian smoothing method        
tibiaV = patchSmooth(tibiaF,tibiaV,[],smoothPar);
fibulaV = patchSmooth(fibulaF,fibulaV,[],smoothPar);

%Remesh using ggremesh to reduce complexity
[tibiaF,tibiaV] = ggremesh(tibiaF, tibiaV, optionStruct_tib);

%Check for matching point number and display if not
if length(tibiaV) ~= optionStruct_tib.nb_pts
    error('Number of tibia vertices does not match requested for case %s\n',caseNo);
end

%Remesh using ggremesh to reduce complexity
[fibulaF,fibulaV] = ggremesh(fibulaF,fibulaV, optionStruct_fib);

%Check for matching point number and display if not
if length(fibulaV) ~= optionStruct_fib.nb_pts
    %Print out summary this time
    error('Number of fibula vertices does not match requested for case %s\n',caseNo);
end

%Check for holes in surfaces
%This shouldn't occur with the remeshing, but check for safety...
if ~isempty(patchBoundary(tibiaF, tibiaV))
    error('Holes in tibia mesh')
end
if ~isempty(patchBoundary(fibulaF, fibulaV))
    error('Holes in fibula mesh')
end

% % % %Visualise imported surfaces
% % % cFigure; hold on
% % % gpatch(tibiaF,tibiaV,'gw','k')
% % % gpatch(fibulaF,fibulaV,'bw','k')
% % % axisGeom; camlight headlight

%Rotate the landmarks to be in the tibial coordinate system like the shape model
[sampleLandmarks, tibiaV, fibulaV] = ...
    alignLandmarksSurfaces(sampleLandmarks, tibiaV, fibulaV);

% % % %Visualise rotated surfaces and landmarks
% % % cFigure; hold on;
% % % %Tibia
% % % plotV(tibiaV, 'g.')
% % % %Fibula
% % % plotV(fibulaV, 'b.')
% % % %Image landmarks
% % % currLandmarks = fieldnames(imageLandmarks);
% % % for ff = 1:length(currLandmarks)
% % %     plotV(sampleLandmarks.(currLandmarks{ff}), 'r.', 'MarkerSize', 25)
% % % end
% % % %Axis parameters
% % % axisGeom;

%% Optimise the shape model parameters to fit the dataset

% % % pcScores = zeros(tibiaFibulaShapeModel.retainPCs,1);

%Create index labels for landmarks on the shape model
%Set some current points for each landmark on shape model
shapeModelLandmarks.LM = [4.497, -8.012, 32.14];
shapeModelLandmarks.MM = [-3.459, 6.469, -32.03];
shapeModelLandmarks.FH = [-1.066, 357.5, 48.19];
shapeModelLandmarks.TT = [30.9, 352.9, -6.18];
shapeModelLandmarks.MC = [-27.9, 383.3, -25.66];
%Identify closest points to those provided to identify index
shapeModelLandmarkNames = fieldnames(shapeModelLandmarks);
for landmarkNo = 1:length(shapeModelLandmarkNames)
    %Calculate distance between points
    ptDist = distancePoints3d(tibiaFibulaShapeModel.meanPoints, ...
        ones(length(tibiaFibulaShapeModel.meanPoints),1) * shapeModelLandmarks.(shapeModelLandmarkNames{landmarkNo}));
    %Find index of minimum and allocate to variable
    shapeModelLandmarkInds.(shapeModelLandmarkNames{landmarkNo}) = find(min(ptDist) == ptDist);
end

[sumError] = calcLandmarkError(pcScores, tibiaFibulaShapeModel, sampleLandmarks, shapeModelLandmarkInds);


% % % %Visualise original vs. reconstructed
% % % cFigure; hold on;
% % % gpatch(tibiaFibulaShapeModel.F, tibiaFibulaShapeModel.meanPoints, 'rw', 'k', 1);
% % % axisGeom;



%% OLD EXAMPLE OF LAB EXPERIMENTAL DATA BELOW...










%% Load and prepare experimental data

%Set participant ID
pID = 'P01';

%Set-up the btk acquisition
c3dAcq = btkReadAcquisition(['samples\',pID,'\',pID,'-static.c3d']);

%Get the markers from the c3d file
c3dMarkers = btkGetMarkers(c3dAcq);

%Extract the desired markers - RLKN, RMKN, RLM, RMM
%Take the average and rename to tibia landmarks system
%(i.e. LC, MC, LM, MM)
sampleLandmarks.LC = mean(c3dMarkers.RLKN);
sampleLandmarks.MC = mean(c3dMarkers.RMKN);
sampleLandmarks.LM = mean(c3dMarkers.RLM);
sampleLandmarks.MM = mean(c3dMarkers.RMM);

%Create new landmarks needed to align tibia to global
sampleLandmarks.IM = midPoint3d(sampleLandmarks.LM,sampleLandmarks.MM);
sampleLandmarks.IC = midPoint3d(sampleLandmarks.LC,sampleLandmarks.MC);

%Load participants tibia and fibula, and associated extracted landmarks
%Tibia
[tibSTLstruct] = import_STL(['samples\',pID,'\',pID,'-tibia.stl']);
tibiaF = tibSTLstruct.solidFaces{1}; %Faces
tibiaV = tibSTLstruct.solidVertices{1}; %Vertices
[tibiaF,tibiaV] = mergeVertices(tibiaF,tibiaV);
%Fibula
[fibSTLstruct] = import_STL(['samples\',pID,'\',pID,'-fibula.stl']);
fibulaF = fibSTLstruct.solidFaces{1}; %Faces
fibulaV = fibSTLstruct.solidVertices{1}; %Vertices
[fibulaF,fibulaV] = mergeVertices(fibulaF,fibulaV);
%Landmarks
landmarksCSV = readtable(['samples\',pID,'\',pID,'-landmarks.csv']);
imageLandmarkNames = landmarksCSV.landmark;
for landmarkNo = 1:length(imageLandmarkNames)
    sampleLandmarks.(imageLandmarkNames{landmarkNo}) = ...
        [landmarksCSV.X(find(strcmp(landmarksCSV.landmark,imageLandmarkNames{landmarkNo}) == 1)), ...
        landmarksCSV.Y(find(strcmp(landmarksCSV.landmark,imageLandmarkNames{landmarkNo}) == 1)), ...
        landmarksCSV.Z(find(strcmp(landmarksCSV.landmark,imageLandmarkNames{landmarkNo}) == 1))];
end

%Create new landmarks needed to align tibia to global
sampleLandmarks.IM = midPoint3d(sampleLandmarks.LM,sampleLandmarks.MM);
sampleLandmarks.IC = midPoint3d(sampleLandmarks.LC,sampleLandmarks.MC);

%Remesh surfaces

%Run initial smoothing on surfaces
smoothPar.n = 2; %2 iterations
smoothPar.Method = 'LAP'; %Laplacian smoothing method        
tibiaV = patchSmooth(tibiaF,tibiaV,[],smoothPar);
fibulaV = patchSmooth(fibulaF,fibulaV,[],smoothPar);

%Remesh using ggremesh to reduce complexity
[tibiaF,tibiaV] = ggremesh(tibiaF, tibiaV, optionStruct_tib);

%Check for matching point number and display if not
if length(tibiaV) ~= optionStruct_tib.nb_pts
    error('Number of tibia vertices does not match requested for case %s\n',caseNo);
end

%Remesh using ggremesh to reduce complexity
[fibulaF,fibulaV] = ggremesh(fibulaF,fibulaV, optionStruct_fib);

%Check for matching point number and display if not
if length(fibulaV) ~= optionStruct_fib.nb_pts
    %Print out summary this time
    error('Number of fibula vertices does not match requested for case %s\n',caseNo);
end

%Check for holes in surfaces
%This shouldn't occur with the remeshing, but check for safety...
if ~isempty(patchBoundary(tibiaF, tibiaV))
    error('Holes in tibia mesh')
end
if ~isempty(patchBoundary(fibulaF, fibulaV))
    error('Holes in fibula mesh')
end

% % % %Visualise imported surfaces
% % % cFigure; hold on
% % % gpatch(tibiaF,tibiaV,'gw','k')
% % % gpatch(fibulaF,fibulaV,'bw','k')
% % % axisGeom; camlight headlight

%Rotate the landmarks to be in the tibial coordinate system like the shape model
[imageLandmarks, imageLandmarks, tibiaV, fibulaV] = ...
    alignLandmarksSurfaces(imageLandmarks, imageLandmarks, tibiaV, fibulaV);






































