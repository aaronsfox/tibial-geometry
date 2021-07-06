function createShapeModel()

%%%%% TODO: add inputs
%%%%% - processSurfaces = True/False (load and process vs. use previous)
%%%%% - nonRigidReg = True/False (do process or use previously registered)
%%%%% - rigidReg = True/False (do process or use previously registered)

processSurfaces = true;
nonRigidReg = true;
rigidReg = true;

%% 

%%%%% TODO: add notes...

% Function to take original segmentations, align to global coordinate,
% non-rigidly register to a master target mesh, rigidly register to master
% target mesh, create shape model

%% Set-up

%Set home directory
homeDir = cd;

%Add supplementary code folders to path
addpath([pwd,'\ShapeModelBuilder']);

%Navigate to segmentation directory
cd('..\Segmentation\');

%Grab the case names
f = dir();
for ff = 3:length(f)
    caseID{ff-2} = f(ff).name;
end
clear ff

%Set options for remeshing
% % % optionStruct.pointSpacing = 2; %Set desired point spacing
optionStruct.nb_pts = 10000; %Set desired number of points
optionStruct.disp_on = 0; % Turn off command window text display

%%%%% TODO: matching the number of points seems to work best for the shape
%%%%% model registration processes - but can't achieve the desired number
%%%%% with the high density (i.e. high number of nodes) surface meshes

%% Load and remesh surfaces

%Check whether to load and process surfaces or take from pre-processed directory
if processSurfaces

    %Loop through cases
    for ii = 1:length(caseID)

        %Navigate to case ID
        cd(caseID{ii});

        %Get case number
        caseNo = strsplit(caseID{ii},'-');
        caseNo = caseNo{2};

        %Load the tibia
        [tibiaSTLstruct] = import_STL([caseNo,'-tibia-cortical-remesh.stl']);
        tibiaF = tibiaSTLstruct.solidFaces{1}; %Faces
        tibiaV = tibiaSTLstruct.solidVertices{1}; %Vertices
        [tibiaF,tibiaV] = mergeVertices(tibiaF,tibiaV);

        %Try to remesh and continue with case if possible
        try

            %Remesh using ggremesh to reduce complexity
            [tibiaF,tibiaV] = ggremesh(tibiaF, tibiaV, optionStruct);

            %Check for desired number of points
            if length(tibiaV) ~= optionStruct.nb_pts
                
% % %                 %Use point cloud downsample combined with alphashape to
% % %                 %reconstruct the mesh with the desired points
% % % 
% % %                 
% % %                 %Get percentage to downsample by based on n points
% % %                 per = optionStruct.nb_pts/length(tibiaV);
% % %                 
% % %                 %Create the downsampled point cloud from the remeshed surface
% % %                 ptCloud = pcdownsample(pointCloud(tibiaV),'random',per);
% % %                 
% % %                 %Create the shape from the new point cloud
% % %                 aShape = alphaShape(ptCloud.Location);
% % %                 
% % %                 %Remesh to the desired points with the new points and faces
% % %                 [F2,V2] = ggremesh(F, V, optionStruct);
% % %                 
% % %                 [tibiaF2,tibiaV2] = ggremesh(tibiaF, tibiaV, optionStruct2);
% % %                 [tibiaF3,tibiaV3] = ggremesh(tibiaF2, tibiaV2, optionStruct);
% % %                 
% % %                 [Fb,Vb]=triRemeshLabel(aShape.boundaryFacets, ptCloud.Location,3);
% % %                 [Fbg,Vbg] = ggremesh(Fb, Vb, optionStruct);
                
                
                errorOut = 'incorrectPts';
                error('Desired number of points not achieved during remeshing.');
            end

            %Check if tibial landmarks are available and align segment if so
            if isfile('LC.txt') && isfile('LM.txt') && isfile('MC.txt') && isfile('MM.txt')

                %Import tibial landmarks
                tibPoints = [{'LC'},{'LM'},{'MC'},{'MM'}];
                for pp = 1:length(tibPoints)
                    tree = xml_read([tibPoints{pp},'.txt']);
                    landmarks.(char(tree.Point.Name)) = tree.Point.Coordinate;% / 1000;
                    clear tree 
                end
                clear pp

                %Create new landmarks
                landmarks.IM = midPoint3d(landmarks.LM,landmarks.MM);
                landmarks.IC = midPoint3d(landmarks.LC,landmarks.MC);

                %Create and align tibial coordinate system

                %Create the tibial planes
                %Frontal
                planes.frontal = createPlane(landmarks.IM,...
                    landmarks.LC,landmarks.MC);
                %Torsional
                planes.torsional = createPlane(landmarks.IC,...
                    landmarks.MM,landmarks.LM);

                %Create transform to get tibia aligned to the global plane
                globalTransform = createBasisTransform3d('global',planes.torsional);

                %Transform surfaces
                for pp = 1:length(tibiaV)
                    tibiaV(pp,:) = transformPoint3d(tibiaV(pp,:),globalTransform);    
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
                for pp = 1:length(tibiaV)
                    tibiaV(pp,:) = transformPoint3d(tibiaV(pp,:),rotZ);    
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
                for pp = 1:length(tibiaV)
                    tibiaV(pp,:) = transformPoint3d(tibiaV(pp,:),transMatrix);    
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
                for pp = 1:length(tibiaV)
                    tibiaV(pp,:) = transformPoint3d(tibiaV(pp,:),createRotationOz(deg2rad(-90)));    
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
                for pp = 1:length(tibiaV)
                    tibiaV(pp,:) = transformPoint3d(tibiaV(pp,:),createRotationOy(deg2rad(-90)));    
                end
                clear pp

                %Transform landmarks
                currLandmarks = fieldnames(landmarks);
                for ff = 1:length(currLandmarks)
                    landmarks.(currLandmarks{ff}) = transformPoint3d(landmarks.(currLandmarks{ff}),createRotationOy(deg2rad(-90)));
                end
                clear ff

% % %                 %Visualise
% % %                 h = cFigure; hold on
% % %                 gpatch(tibiaF,tibiaV,'gw','k');
% % %                 axisGeom; camlight headlight
% % %                 title('Imported and Aligned Tibia');
% % %     
% % %                 %Close figure
% % %                 close all

            else

                %No landmarks present
                fprintf('No landmarks present for case %s. Leaving surface in original position.',caseNo);

            end

            %Export the STL

            %Create the stl structure    
            stlStruct.solidNames = {['tibia-',caseNo]}; %names of parts
            stlStruct.solidVertices = {tibiaV}; %Vertices
            stlStruct.solidFaces = {tibiaF}; %Faces
            stlStruct.solidNormals={[]};

            %Export the stl
            export_STL_txt(['..\..\ShapeModel\segmentations\',caseNo,'-tibia-cortical.stl'], stlStruct);

            %Store in structure
            segmentations.(['case_',caseNo]).F = tibiaF;
            segmentations.(['case_',caseNo]).V = tibiaV;

            %Disp progress
            disp(['Case ',caseNo,' successfully remeshed and exported.']);

        catch

            %Display catch issue
            if strcmp(errorOut,'incorrectPts')
                disp(['Unable to remesh case ',caseNo,' to desired number of points. Skipped...'])
            end

        end

        %Clear variables for re-loop
        clearvars -except caseID homeDir ii optionStruct segmentations processSurfaces nonRigidReg rigidReg

        %Navigate back to segmentation directory
        cd('..');

    end
    clear ii

else
    
    %%%%% TODO: add steps for loading from segmentations directory
    
end

%% Non-rigid registration

%%%%% TODO: this does scaling though???

%Check whether to perform non rigid registration or take from pre-processed directory
if nonRigidReg
    
    %Settings for non rigid registration
    nIters = 10;
    preAlligned = 0;
    
    %Get the case names from the segmentations structure
    caseNames = fieldnames(segmentations);
    
    %Take the first case as the target and place in structure
    fittedMeshes.(caseNames{1}).F = segmentations.(caseNames{1}).F;
    fittedMeshes.(caseNames{1}).V = segmentations.(caseNames{1}).V;
    
    %Loop through remaining cases and register
    for ii = 2:length(caseNames)
        
        %Use ICP algorithm for non rigid registration
        registered = nonrigidICPv1(fittedMeshes.(caseNames{1}).V,....
            segmentations.(caseNames{ii}).V,...
            fittedMeshes.(caseNames{1}).F,...
            segmentations.(caseNames{ii}).F,...
            nIters, preAlligned);
        
        %Close created plot
        close all
        
% % %         %Visualise
% % %         cFigure;
% % %         gpatch(segmentations.(caseNames{ii}).F,registered);
% % %         axisGeom;
        
        %Store in fitted mesh structure
        fittedMeshes.(caseNames{ii}).F = segmentations.(caseNames{ii}).F;
        fittedMeshes.(caseNames{ii}).V = registered;
        
        %Cleanup
        clear registered
        
    end
    clear ii
    
    %Loop through cases and save as surface meshes
    for ii = 2:length(caseNames)
        
        %Get case number
        caseNo = strsplit(caseNames{ii},'_');
        caseNo = caseNo{2};
        
        %Export the STL

        %Create the stl structure    
        stlStruct.solidNames = {['tibia-',caseNo]}; %names of parts
        stlStruct.solidVertices = {fittedMeshes.(caseNames{ii}).V}; %Vertices
        stlStruct.solidFaces = {fittedMeshes.(caseNames{ii}).F}; %Faces
        stlStruct.solidNormals={[]};

        %Export the stl
        export_STL_txt(['..\ShapeModel\fittedMeshes\',caseNo,'-tibia-cortical_nonRigidReg.stl'], stlStruct);
        
    end
    clear ii
    
else
    
    %%%%% TODO: add process for loading non-rigid registered meshes
    
end

%% Rigid registration

%% Shape model

%%%% TODO: change to after rigid???

%Place data in XYZ structures
for ii = 1:length(caseNames)
    xData(:,ii) = fittedMeshes.(caseNames{ii}).V(:,1);
    yData(:,ii) = fittedMeshes.(caseNames{ii}).V(:,2);
    zData(:,ii) = fittedMeshes.(caseNames{ii}).V(:,3);
end
clear ii

%Create shape model
[ssmV,Eval,Evec,MEAN,PCcum,Modes] = SSMbuilder(xData,yData,zData);




%%


%% ----- End of createShapeModel.m ----- %%