


%%%%% TODO: consider denser meshes --- particularly for trabecular???
%%%%% TODO: consider less dense mesh so points don't get mixed up???
%%%%% TODO: consider smoothing post registration???


%%%%% TODO: not sure tibia only shape model works exactly right when the
%%%%% rigid alignment is done on the tib-fib complex
    %%%%% may be the case that we need paired trabecular models with each
    %%%%% so that they align --- could test with tib-fib model to see if
    %%%%% that shape model describes better
    
    
%%%%% TODO: starting to think that a reduced mesh density might actually
%%%%% help out with the shape modelling...

%%%%% TODO: something not entirely quite right with the shape model created
%%%%% --- perhaps it's due to the lack of non-rigid registration and using
%%%%% the original points with that (i.e. maybe we do need the CDP
%%%%% non-rigid registration first and then identify the corresponding
%%%%% points...?)


%% Set-up

%Set home directory
homeDir = cd;

%Add supplementary code path
addpath(genpath([pwd,'\Supplementary']));

%Navigate to segmentation directory
cd('..\Segmentation\');

%Grab the case names
f = dir();
for ff = 3:length(f)
    caseID{ff-2} = f(ff).name;
end

%Set options for remeshing
% % % optionStruct_tib.pointSpacing = 2.5; %Set desired point spacing -- another approach for remeshing
% % % optionStruct_tib.nb_pts = 5000; %Set desired number of points
optionStruct_tib.nb_pts = 4000; %Set desired number of points
optionStruct_tib.disp_on = 0; % Turn off command window text display
% % % optionStruct_fib.pointSpacing = 2.5; %Set desired point spacing -- another approach for remeshing
% % % optionStruct_fib.nb_pts = 2500; %Set desired number of points
optionStruct_fib.nb_pts = 2000; %Set desired number of points
optionStruct_fib.disp_on = 0; % Turn off command window text display

%Set options for CPD algorithm
optCPD.method = 'nonrigid'; % use nonrigid registration
optCPD.beta = 2;            % the width of Gaussian kernel (smoothness)
optCPD.lambda = 3;          % regularization weight
optCPD.viz = 0;             % don't visualise
optCPD.outliers = 0;        % don't account for outliers
optCPD.fgt = 0;             % do not use FGT (default)
optCPD.normalize = 1;       % normalize to unit variance and zero mean before registering (default)
optCPD.corresp = 0;         % compute correspondence vector at the end of registration (not being estimated by default)
optCPD.max_it = 100;        % max number of iterations
optCPD.tol = 1e-4;          % tolerance

%Set options for visualising each remeshing
visualiseRemesh = false; %change to true to view each iteration

%% Remesh and align surfaces

%All surface meshes need to be aligned to a single target surface. For this
%we use a case that is closest to the mean surface in line with existing
% work (i.e. Bruce et al. 2021, Computer Methods Biomech Biomed Eng, doi:
%10.1080/10255842.2021.1985111). Some testing revealed that the best use
%case for this was case-147211, hence we ensure this is processed first and
%all subsequent meshes are aligned to this

%Set reference case number
refCase = '147211';

%Reorder case list to make target mesh first
caseID = [caseID(find(contains(caseID, refCase))), ...
    caseID(1:find(contains(caseID, refCase))-1), ...
    caseID(find(contains(caseID, refCase))+1:end)];

%Remesh and align all cases

%Create waitbar to monitor progress
wbar = waitbar(0,'Remeshing and aligning surfaces....');

%Loop through cases
%%%%% TODO: does MATLAB slow down as the dataset fills up???
tic %%%%%%%% TODO: fix...
for caseInd = 1:length(caseID)
    
    %Navigate to case ID
    cd(caseID{caseInd});
    
    %Get case number
    caseNo = strsplit(caseID{caseInd},'-');
    caseNo = caseNo{2};
    
    %Create the check for no fibula in the case of case-134065
    if strcmp(caseNo, '134065')
        noFib = true;
    else
        noFib = false;
    end
    
    %Load the tibia - cortical
    %Uses GIBBON functionality
    [tibiaSTLstruct] = import_STL([caseNo,'-tibia-cortical-remesh.stl']);
    tibiaF = tibiaSTLstruct.solidFaces{1}; %Faces
    tibiaV = tibiaSTLstruct.solidVertices{1}; %Vertices
    [tibiaF,tibiaV] = mergeVertices(tibiaF,tibiaV);
    
    %Load the tibia - trabecular
    %Uses GIBBON functionality
    [trabSTLstruct] = import_STL([caseNo,'-tibia-trabecular.stl']);
    trabF = trabSTLstruct.solidFaces{1}; %Faces
    trabV = trabSTLstruct.solidVertices{1}; %Vertices
    [trabF,trabV] = mergeVertices(trabF,trabV);
    
    %Load the fibula
    if ~noFib
        %Uses GIBBON functionality
        [fibulaSTLstruct] = import_STL([caseNo,'-fibula.stl']);
        fibulaF = fibulaSTLstruct.solidFaces{1}; %Faces
        fibulaV = fibulaSTLstruct.solidVertices{1}; %Vertices
        [fibulaF,fibulaV] = mergeVertices(fibulaF,fibulaV);
    end

    %Remesh surfaces
        
    %Run initial smoothing on surfaces
    smoothPar.n = 2; %2 iterations
    smoothPar.Method = 'LAP'; %Laplacian smoothing method        
    tibiaV = patchSmooth(tibiaF,tibiaV,[],smoothPar);
    trabV = patchSmooth(trabF,trabV,[],smoothPar);
    if ~noFib
        fibulaV = patchSmooth(fibulaF,fibulaV,[],smoothPar);
    end

    %Remesh using ggremesh to reduce complexity
    [tibiaF,tibiaV] = ggremesh(tibiaF, tibiaV, optionStruct_tib);

    %Check for matching point number and display if not
    if length(tibiaV) ~= optionStruct_tib.nb_pts
        error('Number of tibia vertices does not match requested for case %s\n',caseNo);
    end

    %Remesh using ggremesh to reduce complexity
    [trabF,trabV] = ggremesh(trabF,trabV, optionStruct_tib);

    %Check for matching point number and display if not
    if length(trabV) ~= optionStruct_tib.nb_pts
        %First time here at least retry the remesh again
        [trabF,trabV] = ggremesh(trabF,trabV, optionStruct_tib);
        %Check again
        if length(trabV) ~= optionStruct_tib.nb_pts
            %Print out summary this time
            error('Number of trabecular vertices does not match requested for case %s\n',caseNo);
        end
    end

    %Remesh using ggremesh to reduce complexity
    if ~noFib
        [fibulaF,fibulaV] = ggremesh(fibulaF,fibulaV, optionStruct_fib);

        %Check for matching point number and display if not
        if length(fibulaV) ~= optionStruct_fib.nb_pts
            %First time here at least retry the remesh again
            [fibulaF,fibulaV] = ggremesh(fibulaF,fibulaV, optionStruct_fib);
            %Check again
            if length(fibulaV) ~= optionStruct_fib.nb_pts
                %Print out summary this time
                error('Number of fibula vertices does not match requested for case %s\n',caseNo);
            end
        end
    end

    %Check for holes in surfaces
    %This shouldn't occur with the remeshing, but check for safety...
    if ~isempty(patchBoundary(tibiaF, tibiaV))
        error('Holes in tibia mesh for %s\n...',caseNo)
    end
    if ~isempty(patchBoundary(trabF, trabV))
        error('Holes in trabecular mesh for %s\n...',caseNo)
    end
    if ~noFib
        if ~isempty(patchBoundary(fibulaF, fibulaV))
            error('Holes in fibula mesh for %s\n...',caseNo)
        end
    end

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
    if ~noFib
        [tibiaF, tibiaV, fibulaF, fibulaV, trabF, trabV] = ...
            rotateSurfacesGlobal(landmarks, noFib, tibiaF, tibiaV, fibulaF, fibulaV, trabF, trabV);
    elseif noFib
        [tibiaF, tibiaV, ~, ~, trabF, trabV] = ...
            rotateSurfacesGlobal(landmarks, noFib, tibiaF, tibiaV, [], [], trabF, trabV);
    end

    %Visualise
    if visualiseRemesh
        %Create figure
        h = cFigure; hold on
        %Plot surfaces
        gpatch(tibiaF,tibiaV,'gw','k', 0.3); %low alpha to see trabecular 
        gpatch(trabF,trabV,'rw','k');
        if ~noFib
            gpatch(fibulaF,fibulaV,'bw','k');
        end
        axisGeom; camlight headlight
        title('Imported and Aligned Tibia-Fibula');
        %Pause
        pause
        %Close figure
        close all
    end

    %If not the target mesh, align the surfaces
    if caseInd > 1
        
        %Tibia
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['APPLYING CPD ALGORITHM TO TIBIA OF CASE-',caseNo,'...']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Perform non-rigid registration using CPD to align surfaces
        [cpdTformTib] = cpd_register(data.tibia.V(:,:,1), tibiaV, optCPD);
        
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %         disp(['APPLYING ICP ALGORITHM TO TIBIA OF CASE-',caseNo,'...']);
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %         
% % %         %%%% TODO: could increase to 20 iterations?????
% % %         
% % %         %Perform non-rigid registration using the ICP algorithm to assist in
% % %         %identifying corresponding points.
% % % % % %         [registeredPtsTib] = nonrigidICPv1(data.tibia.V(:,:,1), movingRegTib_cpd.Location, ...
% % % % % %             data.tibia.F(:,:,1), tibiaF, 10, 1);
% % %         [registeredPtsTib] = nonrigidICPv1(data.tibia.V(:,:,1), tibiaV, ...
% % %             data.tibia.F(:,:,1), tibiaF, 10, 0);

        %%%%% TODO: what about knnsearch finding the same points as nearest
        %%%%% neighbours? need to include distance and find closest for
        %%%%% each point without including same? Reducing mesh density
        %%%%% might help with this...?

        %Identify the matching points in the target mesh against the
        %registration to identify the corresponding point indices
        regSortIdxTib = knnsearch(cpdTformTib.Y, data.tibia.V(:,:,1));
        
        %Sort the registered points so they align with the reference mesh

        %%%% This seems appropriate --- but need to consider that fibula
        %%%% and trabecular are being rigidly aligned according to the
        %%%% tibia --- appropriate???
        
% % %         tibiaV_sorted = movingRegTib_cpd.Location(regSortIdxTib,:);
        tibiaV_sorted = tibiaV(regSortIdxTib,:);
        
        %Visualise to confirm that points now correspond. This can be done
        %by using the target mesh faces against the registered moving
        %points. The surface meshes should essentially come out looking the
        %same (i.e. overlapped) despite using different faces
% % %         cFigure; hold on
% % %         gpatch(data.tibia.F(:,:,1), data.tibia.V(:,:,1), 'gw', 'k', 0.3);
% % %         gpatch(tibiaF, tibiaV, 'rw', 'k', 0.3);
% % %         gpatch(data.tibia.F(:,:,1), tibiaV_sorted, 'bw', 'k');
% % %         axisGeom;
        
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %         disp(['APPLYING PROCRUSTES TO TIBIA OF CASE-',caseNo,'...']);
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %         
% % %         %Perform rigid transformation using procrustes analysis
% % %         [tformTib_proc, movingRegTib_proc] = procrustes(data.tibia.V(:,:,1), ...
% % %             tibiaV_sorted, 'scaling', false);
        
        %Allocate to original variables
        tibiaF = data.tibia.F(:,:,1); %use target mesh faces
        tibiaV = tibiaV_sorted; %use sorted registered points
        
        %Trabecular
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['APPLYING CPD ALGORITHM TO TRABECULAR OF CASE-',caseNo,'...']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Perform non-rigid registration using CPD to align surfaces
        [cpdTformTrab] = cpd_register(data.trab.V(:,:,1), trabV, optCPD);
        
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %         disp(['APPLYING ICP ALGORITHM TO TRABECULAR OF CASE-',caseNo,'...']);
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %         
% % %         %Perform non-rigid registration using the ICP algorithm to assist in
% % %         %identifying corresponding points.
% % %         [registeredPtsTrab] = nonrigidICPv1(data.trab.V(:,:,1), trabV, ...
% % %             data.trab.F(:,:,1), trabF, 10, 0);

        %Identify the matching points in the target mesh against the
        %registration to identify the corresponding point indices
        regSortIdxTrab = knnsearch(cpdTformTrab.Y, data.trab.V(:,:,1));
        
        %Sort the registered points so they align with the shape model mean points
        trabV_sorted = trabV(regSortIdxTrab,:);
        
        %Visualise to confirm that points now correspond. This can be done
        %by using the target mesh faces against the registered moving
        %points. The surface meshes should essentially come out looking the
        %same (i.e. overlapped) despite using different faces
% % %         cFigure; hold on
% % %         gpatch(data.trab.F(:,:,1), data.trab.V(:,:,1), 'gw', 'k', 0.3);
% % %         gpatch(trabF, trabV, 'rw', 'k', 0.3);
% % %         gpatch(data.trab.F(:,:,1), trabV_sorted, 'bw', 'k');
% % %         axisGeom;
        
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %         disp(['APPLYING PROCRUSTES TO TRABECULAR OF CASE-',caseNo,'...']);
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %         
% % %         %Perform rigid transformation using procrustes analysis
% % %         [tformTrab_proc, movingRegTrab_proc] = procrustes(data.trab.V(:,:,1), ...
% % %             trabV_sorted, 'scaling', false);
        
        %Allocate to original variables
        trabF = data.trab.F(:,:,1); %use target mesh faces
        trabV = trabV_sorted; %use sorted registered points
        
        %Fibula

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['APPLYING CPD ALGORITHM TO FIBULA OF CASE-',caseNo,'...']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %Perform non-rigid registration using CPD to align surfaces
        [cpdTformFib] = cpd_register(data.fibula.V(:,:,1), fibulaV, optCPD);

% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %         disp(['APPLYING ICP ALGORITHM TO FIBULA OF CASE-',caseNo,'...']);
% % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % %         %Perform non-rigid registration using the ICP algorithm to assist in
% % %         %identifying corresponding points.
% % %         [registeredPtsFib] = nonrigidICPv1(data.fibula.V(:,:,1), fibulaV, ...
% % %             data.fibula.F(:,:,1), fibulaF, 10, 0);

        %Identify the matching points in the target mesh against the
        %registration to identify the corresponding point indices
        regSortIdxFib = knnsearch(cpdTformFib.Y, data.fibula.V(:,:,1));

        %Sort the registered points so they align with the shape model mean points
        fibulaV_sorted = fibulaV(regSortIdxFib,:);

        %Visualise to confirm that points now correspond. This can be done
        %by using the target mesh faces against the registered moving
        %points. The surface meshes should essentially come out looking the
        %same (i.e. overlapped) despite using different faces
% % %         cFigure; hold on
% % %         gpatch(data.fibula.F(:,:,1), data.fibula.V(:,:,1), 'gw', 'k', 0.3);
% % %         gpatch(fibulaF, fibulaV, 'rw', 'k', 0.3);
% % %         gpatch(data.fibula.F(:,:,1), fibulaV_sorted, 'bw', 'k');
% % %         axisGeom;

% % %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %             disp(['APPLYING PROCRUSTES TO TIBIA OF CASE-',caseNo,'...']);
% % %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % %             %Perform rigid transformation using procrustes analysis
% % %             [tformFib_proc, movingRegFib_proc] = procrustes(data.fibula.V(:,:,1), ...
% % %                 fibulaV_sorted, 'scaling', false);

        %Allocate to original variables
        fibulaF = data.fibula.F(:,:,1); %use target mesh faces
        fibulaV = fibulaV_sorted; %use sorted registered points

        %Merge the tibia and fibula vertices for rigid alignment
        base_tibia_fibula = [data.tibia.V(:,:,1);data.fibula.V(:,:,1)];
        tibia_fibulaV = [tibiaV;fibulaV];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['APPLYING PROCRUSTES TO SURFACES OF CASE-',caseNo,'...']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%% TODO: what are the implications for the different shape
        %%%%% models for applying transformation based on tibia vs.
        %%%%% tibia-fibula...?

        %Perform rigid transformation using procrustes analysis
        [~, regV, tform] = procrustes(base_tibia_fibula, ...
            tibia_fibulaV, 'scaling', false);
% % %         [~, regV, tform] = procrustes(data.tibia.V(:,:,1), ...
% % %             tibiaV, 'scaling', false);

        %Split data back to individual surfaces
        tibiaV = regV(1:optionStruct_tib.nb_pts,:);
        fibulaV = regV(optionStruct_tib.nb_pts+1:end,:);
% % %         tibiaV_proc = regV;

        %Apply the rigid transform to the trabecular for alignment
        trabV = trabV * tform.T + mean(tform.c);
% % %         fibulaV_proc = fibulaV * tform.T + mean(tform.c);
        
        %Visualise to confirm if need be
% % %         cFigure; hold on
% % %         %Plot surfaces
% % %         gpatch(tibiaF,tibiaV_proc,'gw','k', 0.3); %low alpha to see trabecular 
% % %         gpatch(trabF,trabV_proc,'rw','k');
% % %         gpatch(fibulaF,fibulaV_proc,'bw','k');
        
        %%%%% TODO: export subplots of registrations to go back and review
        %%%%% visually...
        
    end
    
    %Store data in structure
    
    %Tibia models
    data.tibia.caseID{caseInd} = caseID{caseInd};
    data.tibia.F(:,:,caseInd) = tibiaF;
    data.tibia.V(:,:,caseInd) = tibiaV;
  
    %Trabecular models
    data.trab.caseID{caseInd} = caseID{caseInd};
    data.trab.F(:,:,caseInd) = trabF;
    data.trab.V(:,:,caseInd) = trabV;
    
    %Fibula models
    data.fibula.caseID{caseInd} = caseID{caseInd};
    data.fibula.F(:,:,caseInd) = fibulaF;
    data.fibula.V(:,:,caseInd) = fibulaV;
    
    %Update waitbar
    waitbar(caseInd/length(caseID), wbar, ...
        ['Remeshed and aligned ',num2str(caseInd),' of ',num2str(length(caseID)),' surfaces...']);
    
    %Disp progress
    disp(['Case ',caseNo,' successfully remeshed and aligned.']);
    
    %Clear variables for re-loop
    clearvars -except visualiseRemesh caseID homeDir caseInd optionStruct_tib optionStruct_fib optCPD wbar data 
    
    %Navigate back to segmentation directory
    cd('..');
    
end

%Update and close waitbar
waitbar(1, wbar, 'Remeshed all surfaces. Finishing...');
close(wbar);
duration = toc


%%


%% Create shape models

%% Tibia only shape model

%Loop through case ID's and reshape data for PCA
for caseInd = 1:length(caseID)
    
    %Reshape data
    nodes(caseInd,:) = reshape(data.tibia.V(:,:,caseInd)', optionStruct_tib.nb_pts*3, [])';
    
end

%Run PCA to create shape model
[shapeModel.loadings, shapeModel.score, shapeModel.latent, ~, ...
    shapeModel.varExplained, shapeModel.mean] = pca(nodes);


% % % cFigure; hold on
% % % mnPts = reshape(shapeModel.mean, [3, length(shapeModel.mean)/3])';
% % % plotV(mnPts,'b.');
% % % gpatch(data.tibia.F(:,:,1),mnPts,'rw','k');
% % % axisGeom;

%% Tibia-fibula shape model

%Loop through case ID's and reshape data for PCA
for caseInd = 1:length(caseID)
    
    %Reshape data
    nodes(caseInd,:) = reshape([data.tibia.V(:,:,caseInd); ...
        data.fibula.V(:,:,caseInd)]',...
        (optionStruct_tib.nb_pts + optionStruct_fib.nb_pts)*3, [])';
    
end

%Run PCA to create shape model
[shapeModel.loadings, shapeModel.score, shapeModel.latent, ~, ...
    shapeModel.varExplained, shapeModel.mean] = pca(nodes);

%Create new faces structure to visualise
F = data.tibia.F(:,:,1);
% % % F = [data.tibia.F(:,:,1); data.fibula.F(:,:,1) + length(data.tibia.V)];


% % % cFigure; hold on
% % % mnPts = reshape(shapeModel.mean, [3, length(shapeModel.mean)/3])';
% % % plotV(mnPts,'b.');
% % % gpatch(F,mnPts,'rw','k');
% % % axisGeom;


%% Attempt to reconstruct case surface from principal components
%  Aim here is to test out process of creating 'new' models

%Reconstruct the mean centred data for each participant
reconstructed = shapeModel.score * shapeModel.loadings';

%Test out a specific case
caseInd = 5;
%Extract data
testCase = reconstructed(caseInd,:);
%Add the mean and reshape
reconstructedV = reshape(testCase + shapeModel.mean, [3, length(shapeModel.mean)/3])';

%Visualise to compare
cFigure; hold on;
%Original
gpatch(data.tibia.F(:,:,caseInd), data.tibia.V(:,:,caseInd),'gw','none',0.3);
% % % gpatch(data.fibula.F(:,:,caseInd), data.fibula.V(:,:,caseInd),'gw','none',0.3);
%Reconstructed
% gpatch(data.tibia.F(:,:,caseInd), reconstructedV,'rw','k',0.3);
gpatch(F, reconstructedV,'rw','k');
plotV(reconstructedV,'b.');
%Axes
axisGeom;


%% Create animation of shape change for PC

%Set PC
PC = 5;

%Set number of PCs
reconPCs = length(shapeModel.varExplained);

%Reshape mean
meanPts = reshape(shapeModel.mean, [3, length(shapeModel.mean)/3])';

%Calculate mean and standard deviation values for scores
% % % meanScore = mean(shapeModel.score); %basically zero
sdScore = std(shapeModel.score);

%Create a variable to work through 0.1 increments from -3 to +3 SD
calcSD = linspace(0.1, 3, diff([0,3])*10);

%Calculate the points at each increment
for nSD = 1:length(calcSD)
    
    %Reconstruct the points at the current SD interval

    %Set the simulated scores array with the standard deviation added/subtracted
    %Plus SD
    simScorePlus = zeros(1,reconPCs);
    simScorePlus(PC) = simScorePlus(PC) + (sdScore(PC) * calcSD(nSD));
    %Minus SD
    simScoreMinus = zeros(1,reconPCs);
    simScoreMinus(PC) = simScoreMinus(PC) + (sdScore(PC) * -calcSD(nSD));
    
    %Reconstruct the points, add the mean and reshape to 3D
    plusPts = reshape((simScorePlus * shapeModel.loadings(:,1:reconPCs)') + shapeModel.mean, ...
        [3, length(shapeModel.mean)/3])';
    minusPts = reshape((simScoreMinus * shapeModel.loadings(:,1:reconPCs)') + shapeModel.mean, ...
        [3, length(shapeModel.mean)/3])';
        
    %Set into structure
    posV(:,:,nSD) = plusPts;
    negV(:,:,nSD) = minusPts;
    
    %Calculate distances between points for heat mapping
    pcPlusDist = distancePoints3d(meanPts, plusPts);
    pcMinusDist = distancePoints3d(meanPts, minusPts);

    %Convert distances to face measures
%     pcPlusC = vertexToFaceMeasure(data.tibia.F(:,:,1), pcPlusDist);
%     pcMinusC = vertexToFaceMeasure(data.tibia.F(:,:,1), pcMinusDist);
    pcPlusC = vertexToFaceMeasure(F, pcPlusDist);
    pcMinusC = vertexToFaceMeasure(F, pcMinusDist);

    %Set into structure
    posC(:,nSD) = pcPlusC;
    negC(:,nSD) = pcMinusC; 
    
end

%Create the basic view to store graphics and initiate animation
hf = cFigure;
%Positive SD
subplot(1,2,1);
title(['Plus SD for PC',num2str(PC)]);
% hp_pos = gpatch(data.tibia.F(:,:,1), meanPts, posC(:,1), 'none', 1);
hp_pos = gpatch(F, meanPts, posC(:,1), 'none', 1);
axisGeom; camlight headlight
axis(axisLim([posV;negV])); %set axis limits based on values
axis off
%Negative SD
subplot(1,2,2);
title(['Minus SD for PC',num2str(PC)]);
% hp_neg = gpatch(data.tibia.F(:,:,1), meanPts, negC(:,1), 'none', 1);
hp_neg = gpatch(F(:,:,1), meanPts, negC(:,1), 'none', 1);
axisGeom; camlight headlight
axis(axisLim([posV;negV])); %set axis limits based on values
axis off

%Set up animation
animStruct.Time = calcSD;
%Loop through SD points
for nSD = 1:length(calcSD)
    
    %Set entries into animation structure
    animStruct.Handles{nSD} = [hp_pos hp_pos hp_neg hp_neg]; %Handles of objects to animate
    animStruct.Props{nSD} = {'Vertices', 'CData', 'Vertices', 'CData'}; %Properties of objects to animate
    animStruct.Set{nSD} = {posV(:,:,nSD), posC(:,nSD), negV(:,:,nSD), negC(:,nSD)}; %Property values for to set in order to animate
    
end

%Animate figure
anim8(hf, animStruct);

%%%% TODO: could provide alternative views across subplots in anim8 figure?

%%%% NOTE: for PC1 of tib-fib there's a bit of noise at larger SDs as the
%%%% nodes get squished together or pulled apart --- unsure of how to
%%%% address that??? It occurs across other components too --- it's like a
%%%% distortion --- but that might be as I'm testing it with only 10
%%%% surfaces...
%%%%
%%%% Outside of this, the tibia-fibula model seems to be picking up typical
%%%% shape variation versus just random noise...
%%%%
%%%% Doesn't get much better necessarily with more cases in the shape model
%%%% - but it's perhaps a little more noticeable with the color mapping, as
%%%% coincident points have more color variation







