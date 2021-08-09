%% 

% Function to take original segmentations, align to global coordinate
% system and then output to the starting shape model folder.

%% Set-up

%Set home directory
homeDir = cd;

%Navigate to segmentation directory
cd('..\Segmentation\');

%Grab the case names
f = dir();
for ff = 3:length(f)
    caseID{ff-2} = f(ff).name;
end
clear ff

%Set options for remeshing
% % % optionStruct.pointSpacing = 2.5; %Set desired point spacing
optionStruct_tib.nb_pts = 10000; %Set desired number of points
% % % optionStruct.anisotropy = 1;
optionStruct_tib.disp_on = 0; % Turn off command window text display
optionStruct_fib.nb_pts = 2500; %Set desired number of points
optionStruct_fib.disp_on = 0; % Turn off command window text display

%% Extract and process surfaces

%Loop through cases
% for ii = 1:length(caseID)
caseID = [{'case-102480'};
    {'case-102924'};
    {'case-103559'};
    {'case-103862'};
    {'case-107215'};
    {'case-107813'}];
for ii = 1:length(caseID) %%%%%%%%%% FIX this up after testing!!!!!
    
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
    
    %Load the fibula
    [fibulaSTLstruct] = import_STL([caseNo,'-fibula.stl']);
    fibulaF = fibulaSTLstruct.solidFaces{1}; %Faces
    fibulaV = fibulaSTLstruct.solidVertices{1}; %Vertices
    [fibulaF,fibulaV] = mergeVertices(fibulaF,fibulaV);

    %Try to remesh and continue with case if possible
    try

        %Remesh using ggremesh to reduce complexity
        [tibiaF,tibiaV] = ggremesh(tibiaF, tibiaV, optionStruct_tib);
        
        %Check for matching point number and display if not
        if length(tibiaV) ~= optionStruct_tib.nb_pts
            fprintf('Number of tibia vertices does not match requested for case %s\n',caseNo);
        end
        
        %Remesh using ggremesh to reduce complexity
        [fibulaF,fibulaV] = ggremesh(fibulaF,fibulaV, optionStruct_fib);
        
        %Check for matching point number and display if not
        if length(fibulaV) ~= optionStruct_fib.nb_pts
            fprintf('Number of fibula vertices does not match requested for case %s\n',caseNo);
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
            for pp = 1:length(fibulaV)
                fibulaV(pp,:) = transformPoint3d(fibulaV(pp,:),globalTransform);    
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
            for pp = 1:length(fibulaV)
                fibulaV(pp,:) = transformPoint3d(fibulaV(pp,:),rotZ);    
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
            for pp = 1:length(fibulaV)
                fibulaV(pp,:) = transformPoint3d(fibulaV(pp,:),transMatrix);    
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
            for pp = 1:length(fibulaV)
                fibulaV(pp,:) = transformPoint3d(fibulaV(pp,:),createRotationOz(deg2rad(-90)));    
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
            for pp = 1:length(fibulaV)
                fibulaV(pp,:) = transformPoint3d(fibulaV(pp,:),createRotationOy(deg2rad(-90)));    
            end
            clear pp

            %Transform landmarks
            currLandmarks = fieldnames(landmarks);
            for ff = 1:length(currLandmarks)
                landmarks.(currLandmarks{ff}) = transformPoint3d(landmarks.(currLandmarks{ff}),createRotationOy(deg2rad(-90)));
            end
            clear ff

% % %             %Visualise
% % %             h = cFigure; hold on
% % %             gpatch(tibiaF,tibiaV,'gw','k');
% % %             gpatch(fibulaF,fibulaV,'bw','k');
% % %             axisGeom; camlight headlight
% % %             title('Imported and Aligned Tibia-Fibula');
% % % 
% % %             %Save figure
% % %             saveas(h,'alignedForShapeModel.png')
% % % 
% % %             %Close figure
% % %             close all
            
        else
            
            %No landmarks present
            fprintf('No landmarks present for case %s. Leaving surface in original position.',caseNo);
            
        end

        %Export the STL
        
        %Merge objects into one surface
        mergedSurfaceV = [tibiaV; fibulaV];
        mergedSurfaceF = [tibiaF; fibulaF+length(tibiaV)];

        %Create the stl structure    
        stlStruct.solidNames = {['tibia-fibula-',caseNo]}; %names of parts
        stlStruct.solidVertices = {mergedSurfaceV}; %Vertices
        stlStruct.solidFaces = {mergedSurfaceF}; %Faces
        stlStruct.solidNormals={[]};

        %Export the stl
% % %         export_STL_txt(['..\..\ShapeModel\segmentations\',caseNo,'-tibia-cortical.stl'], stlStruct);
        export_STL_txt(['..\..\ShapeModel_TibFib\segmentations\',caseNo,'-tibia-fibula.stl'], stlStruct);
        
        %Write the case details to a string array for the pca text files
        rbfreg{ii,1} = ['segmentations/',caseNo,'-tibia-fibula.stl'];
        rigidreg{ii,1} = ['fitted_meshes/',caseNo,'-tibia-fibula_rbfreg.stl'];
        pcalist{ii,1} = ['aligned_meshes/',caseNo,'-tibia-fibula_rbfreg_rigidreg.stl'];

        %Disp progress
        disp(['Case ',caseNo,' successfully remeshed and aligned.']);
        
    catch
        
        %Display catch issue
        disp(['Unable to remesh case ',caseNo,'. Skipped...'])
    
    end
    
    %Clear variables for re-loop
    clearvars -except caseID homeDir ii optionStruct_tib optionStruct_fib pcalist rbfreg rigidreg
    
    %Navigate back to segmentation directory
    cd('..');
    
end
clear ii

%Write the case details to the text files needed for the shape model

%Start by removing any old list files
%Remove rbfreg file if present
if isfile('..\ShapeModel_TibFib\rbfreg_list.txt')
    delete('..\ShapeModel_TibFib\rbfreg_list.txt')
end
%Remove rigidreg if present
if isfile('..\ShapeModel_TibFib\rigidreg_list.txt')
    delete('..\ShapeModel_TibFib\rigidreg_list.txt')
end
%Remove pca list if present
if isfile('..\ShapeModel_TibFib\pca_list.txt')
    delete('..\ShapeModel_TibFib\pca_list.txt')
end

%Write the created list sets to file
%rbfreg
fid = fopen('..\ShapeModel_TibFib\rbfreg_list.txt','wt');
for ii = 1:length(rbfreg)
    %Check if current case was successful
    if ~isempty(rbfreg{ii,1})
        %Write to file
        fprintf(fid, rbfreg{ii,1});
        %Check to see if new line is needed
        if ii < length(rbfreg)
            fprintf(fid,'\n');
        end
    end
end
fclose(fid);
%rigidreg
fid = fopen('..\ShapeModel_TibFib\rigidreg_list.txt','wt');
for ii = 1:length(rigidreg)
    %Check if current case was successful
    if ~isempty(rigidreg{ii,1})
        %Write to file
        fprintf(fid, rigidreg{ii,1});
        %Check to see if new line is needed
        if ii < length(rigidreg)
            fprintf(fid,'\n');
        end
    end
end
fclose(fid);
%pca
fid = fopen('..\ShapeModel_TibFib\pca_list.txt','wt');
for ii = 1:length(pcalist)
    %Check if current case was successful
    if ~isempty(pcalist{ii,1})
        %Write to file
        fprintf(fid, pcalist{ii,1});
        %Check to see if new line is needed
        if ii < length(pcalist)
            fprintf(fid,'\n');
        end
    end
end
fclose(fid);

%% Finish up

%Return to home directory
cd(homeDir);

%% ----- End of shapeModelSetUp.m ----- %%