%% Add code notes

%   - Run from current directory

%% Import surfaces

%Tibia
[tibiaSTLstruct] = import_STL('tibia.stl');
tibiaF = tibiaSTLstruct.solidFaces{1}; %Faces
tibiaV = tibiaSTLstruct.solidVertices{1}; %Vertices
[tibiaF,tibiaV] = mergeVertices(tibiaF,tibiaV);

%Fibula
[fibulaSTLstruct] = import_STL('fibula.stl');
fibulaF = fibulaSTLstruct.solidFaces{1}; %Faces
fibulaV = fibulaSTLstruct.solidVertices{1}; %Vertices
[fibulaF,fibulaV] = mergeVertices(fibulaF,fibulaV);

%Joint surface (nodes only)
[jointSurfaceSTLstruct] = import_STL('jointSurface.stl');
surfaceFaces = jointSurfaceSTLstruct.solidFaces{1}; %Faces
surfaceNodes = jointSurfaceSTLstruct.solidVertices{1}; %Vertices
[surfaceFaces,surfaceNodes] = mergeVertices(surfaceFaces,surfaceNodes);

%% Import and create landmarks

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

%% Create and align tibial coordinate system

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
for pp = 1:length(surfaceNodes)
    surfaceNodes(pp,:) = transformPoint3d(surfaceNodes(pp,:),globalTransform);    
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
for pp = 1:length(surfaceNodes)
    surfaceNodes(pp,:) = transformPoint3d(surfaceNodes(pp,:),rotZ);    
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
for pp = 1:length(surfaceNodes)
    surfaceNodes(pp,:) = transformPoint3d(surfaceNodes(pp,:),transMatrix);    
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
for pp = 1:length(surfaceNodes)
    surfaceNodes(pp,:) = transformPoint3d(surfaceNodes(pp,:),createRotationOz(deg2rad(-90)));    
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
for pp = 1:length(surfaceNodes)
    surfaceNodes(pp,:) = transformPoint3d(surfaceNodes(pp,:),createRotationOy(deg2rad(-90)));    
end
clear pp

%Transform landmarks
currLandmarks = fieldnames(landmarks);
for ff = 1:length(currLandmarks)
    landmarks.(currLandmarks{ff}) = transformPoint3d(landmarks.(currLandmarks{ff}),createRotationOy(deg2rad(-90)));
end
clear ff

%% Reduce mesh size for speed

%Remesh surfaces
[tibiaF,tibiaV] = triRemeshLabel(tibiaF,tibiaV,3.0);
[fibulaF,fibulaV] = triRemeshLabel(fibulaF,fibulaV,3.0);

%% Visualise imported and rotated surfaces

%Visualise
cFigure; hold on
gpatch(tibiaF,tibiaV,'gw','k');
gpatch(fibulaF,fibulaV,'bw','k');
axisGeom; camlight headlight
title('Imported and Aligned Surfaces');

%% Cut top of tibia surface for fixed boundary conditions

%Set slicing parameters
snapTolerance = mean(patchEdgeLengths(tibiaF,tibiaV))/100;

%Set normal direction of cutting plane as along the vertical y-axis
n = [0 1 0];

%Set point on the plane as 75% along the length of the tibia
tibiaLength = distancePoints3d(landmarks.IM,landmarks.IC);
p = [0,landmarks.IM(2)+(tibiaLength*0.75),0];

%Slice surface (note 3rd color data output is supressed)
[tibiaFc,tibiaVc,~,logicSide,tibiaEb] = triSurfSlice(tibiaF,tibiaV,[],...
    p,n,snapTolerance);

%Visualise sliced mesh
cFigure;
subplot(1,2,1); hold on;
hp1=gpatch(tibiaFc(~logicSide,:),tibiaVc,'bw','none',1);
hp2=gpatch(tibiaFc(logicSide,:),tibiaVc,'rw','none',1);
legend([hp1 hp2],{'Surface above plane','Surface below plane'})
axisGeom; axis manual; camlight headlight;
colormap gjet;

subplot(1,2,2); hold on;
gpatch(tibiaFc(logicSide,:),tibiaVc,'w','none',1);
gpatch(tibiaFc(~logicSide,:),tibiaVc,'w','none',0.25);
hp1=gpatch(tibiaEb,tibiaVc,'none','b',1,3);
hp2=quiverVec(p,n,50,'k');
legend([hp1 hp2],{'Intersection curve','Plane normal vector'})
axisGeom; axis manual; camlight headlight;

%Extract tibia
tibiaModF = tibiaFc(logicSide==1,:);
[tibiaModF,tibiaModV] = patchCleanUnused(tibiaModF,tibiaVc);

%Fill the hole mad by the cut

%Get the surface boundary
Eb = patchBoundary(tibiaModF,tibiaModV);
indCurve = edgeListToCurve(Eb);
indCurve = indCurve(1:end-1);

%Smooth edge boundary
cparSmooth.n = 5;
cparSmooth.Method = 'HC';
[Eb_smooth] = patchSmooth(Eb,tibiaModV(:,[1 2]),[],cparSmooth);
tibiaModV(indCurve,[1 2]) = Eb_smooth(indCurve,:);
cparSmooth.RigidConstraints = indCurve;
[tibiaModV] = patchSmooth(tibiaModF,tibiaModV,[],cparSmooth);

%Mesh the hole region
pointSpacing = mean(patchEdgeLengths(tibiaModF,tibiaModV));
[fillF,fillV] = regionTriMesh3D({tibiaModV(indCurve,:)},pointSpacing,0,'linear');
%Flip normals if inverted (i.e. -y)
if dot(mean(patchNormal(fillF,fillV)),[0 -1 0]) > 0
    fillF = fliplr(fillF);
end

%Join element sets
[tibiaModF,tibiaModV,tibiaModC] = joinElementSets({tibiaModF,fillF},...
    {tibiaModV,fillV});
[tibiaModF,tibiaModV] = mergeVertices(tibiaModF,tibiaModV);

%% Mesh using tetgen

%Find interior point
innerPoint = getInnerPoint(tibiaModF,tibiaModV);

% % % %Visualise interior point to confirm
% % % cFigure; hold on;
% % % gpatch(tibiaModF,tibiaModV,'w','none',0.5);
% % % plotV(innerPoint,'r.','MarkerSize',25)
% % % axisGeom; camlight headlight;

%Set mesh parameters
tetVolume = tetVolMeanEst(tibiaModF,tibiaModV); %Volume for regular tets
tetGenStruct.stringOpt = '-pq1.2AaY';
tetGenStruct.Faces = tibiaModF;
tetGenStruct.Nodes = tibiaModV;
tetGenStruct.holePoints = [];
tetGenStruct.faceBoundaryMarker = tibiaModC; %Face boundary markers
tetGenStruct.regionPoints = innerPoint; %region points
tetGenStruct.regionA = tetVolume;

%Run tetgen
[meshOutput] = runTetGen(tetGenStruct);

%Access elements, nodes, and boundary faces
tibiaMeshE = meshOutput.elements;
tibiaMeshV = meshOutput.nodes;
tibiaMeshFb = meshOutput.facesBoundary;
tibiaMeshCb = meshOutput.boundaryMarker;
tibiaMeshCE = meshOutput.elementMaterialID;

%Visualise solid mesh
hFig = cFigure; hold on;
optionStruct.hFig = hFig;
meshView(meshOutput,optionStruct);
axisGeom;

%% Identify ankle joint contact surface nodes for prescribed force

%Find the closest nodes on the tibia mesh to the surface nodes
for pp = 1:length(surfaceNodes)
    checkPt = surfaceNodes(pp,:);
    ptDist = distancePoints3d(tibiaModV(:,:),checkPt);
    jointNode(pp,:) = tibiaModV(find(ptDist == min(ptDist)),:);
end
clear pp

%Find index on mesh that match joint node points and corresponding faces
logicJointNodes = ismember(tibiaMeshV,jointNode,'rows');
logicJointFaces = all(logicJointNodes(tibiaMeshFb),2);
bcForcePrescribeList = unique(tibiaMeshFb(logicJointFaces,:));

%Visualise to confirm
cFigure; hold on;
gpatch(tibiaMeshFb,tibiaMeshV,'w','none',1);
plotV(tibiaMeshV(bcForcePrescribeList,:),'r.','markerSize',15)
axisGeom; camlight headlight;

%Note this process could be performed better, as the joint contact surface
%was manually brushed/labelled in 3matic

%% Work out force distribution on joint surface

%Note that this is based on surface normal directions. Forces are assumed
%to only be able to act in a compressive sense on the bone, as per the
%GIBBON example.

%Get nodal normal directions
[~,~,N] = patchNormal(fliplr(tibiaMeshFb),tibiaMeshV);

%Define applied force
%Relatively arbitrary here, but tried to take the axial (Y),
%anterior-posterior (X), and medial-lateral (Z) mean body weight forces
%from Edwards et al. (2010) and convert back to newtons using the mean body
%weight from the same study (70.1kg), using the 4.5m/s running
axialContactF = 13.80;
apContactF = -0.66;
mlContactF = 0.68;
mass = 70.1;
forceTotal = [apContactF * (mass*9.80665), ...
    axialContactF * (mass*9.80665), ...
    mlContactF * (mass*9.80665)];
%%%% Estimate seems somewhat large????

%Set force vectors
FX = [forceTotal(1) 0 0]; %X force vector
FY = [0 forceTotal(2) 0]; %Y force vector
FZ = [0 0 forceTotal(3)]; %Z force vector

%Set force on prescribe nodes
wx = dot(N(bcForcePrescribeList,:),FX(ones(numel(bcForcePrescribeList),1),:),2);
wy = dot(N(bcForcePrescribeList,:),FY(ones(numel(bcForcePrescribeList),1),:),2);
wz = dot(N(bcForcePrescribeList,:),FZ(ones(numel(bcForcePrescribeList),1),:),2);

%Force zero
wx(wx>0) = 0; wy(wy>0) = 0; wz(wz>0) = 0;

%Calculate forces
force_X = forceTotal(1).*ones(numel(bcForcePrescribeList),1).*wx;
force_Y = forceTotal(2).*ones(numel(bcForcePrescribeList),1).*wy;
force_Z = forceTotal(3).*ones(numel(bcForcePrescribeList),1).*wz;

force_X = force_X./sum(force_X(:)); %sum now equal to 1
force_X = force_X.*forceTotal(1); %sum now equal to desired

force_Y = force_Y./sum(force_Y(:)); %sum now equal to 1
force_Y = force_Y.*forceTotal(2); %sum now equal to desired

force_Z = force_Z./sum(force_Z(:)); %sum now equal to 1
force_Z = force_Z.*forceTotal(3); %sum now equal to desired

%Visualise force distributions
%X axis
cFigure;
subplot(1,3,1);hold on;
title('F_x');
gpatch(tibiaMeshFb,tibiaMeshV,'w','none',0.5);
quiverVec([0 0 0],FX,100,'k');
quiverVec(tibiaMeshV(bcForcePrescribeList,:),N(bcForcePrescribeList,:),10,force_X);
axisGeom; camlight headlight;
colormap(gca,gjet(250)); colorbar;
%Y axis
subplot(1,3,2);hold on;
title('F_y');
gpatch(tibiaMeshFb,tibiaMeshV,'w','none',0.5);
quiverVec([0 0 0],FY,100,'k');
quiverVec(tibiaMeshV(bcForcePrescribeList,:),N(bcForcePrescribeList,:),10,force_Y);
axisGeom; camlight headlight;
colormap(gca,gjet(250)); colorbar;
%Z axis
subplot(1,3,3);hold on;
title('F_z');
gpatch(tibiaMeshFb,tibiaMeshV,'w','none',0.5);
quiverVec([0 0 0],FZ,100,'k');
quiverVec(tibiaMeshV(bcForcePrescribeList,:),N(bcForcePrescribeList,:),10,force_Z);
axisGeom; camlight headlight;
colormap(gca,gjet(250)); colorbar;

%%%%% Unsure if forces are being applied in exactly the right directions,
%%%%% so this is something to confirm/check

%% Visualise boundary conditions

%Set support list based on filled upper surface of the tibia
topSupport = tibiaMeshFb(tibiaMeshCb == 2,:);
bcSupportList = unique(topSupport(:));

%Visualise
hFig = cFigure; hold on;
gpatch(tibiaMeshFb,tibiaMeshV,'kw','none',0.25);
hl(1) = plotV(tibiaMeshV(bcSupportList,:),'k.','MarkerSize',25);
hl(2) = plotV(tibiaMeshV(bcForcePrescribeList,:),'r.','MarkerSize',25);
legend(hl,{'BC support','BC force prescribe'});
axisGeom; camlight headlight;

%% Define FEBio input structure

%Settings

%File details
febioFebFileNamePart = 'tibiaModel';
febioFebFileName = fullfile(pwd,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName = fullfile(pwd,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp = [febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_force = [febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force
febioLogFileName_stress = [febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stresses
febioLogFileName_strainEnergy = [febioFebFileNamePart,'_energy_out.txt']; %Log file name for exporting strain energy density

%Material parameters (MPa if spatial units are mm)
%Parameters here are somewhat an average of cortical and trabecular bone.
%In an ideal setting we'd split these up within the mesh
youngsMod = 15000; %Youngs modulus
poisson = 0.3; %Poissons ratio

%FEA control settings
numTimeSteps = 10; %Number of time steps desired
max_refs = 25; %Max reforms
max_ups = 0; %Set to zero to use full-Newton iterations
opt_iter = 6; %Optimum number of iterations
max_retries = 5; %Maximum number of retires
dtmin = (1/numTimeSteps)/100; %Minimum time step size
dtmax = 1/numTimeSteps; %Maximum time step size
runMode = 'internal'; %'external' or 'internal'

%Get a template with default settings
[febio_spec] = febioStructTemplate;

%febio_spec version
febio_spec.ATTR.version = '2.5';

%Module section
febio_spec.Module.ATTR.type = 'solid';

%Control section
febio_spec.Control.analysis.ATTR.type = 'static';
febio_spec.Control.time_steps = numTimeSteps;
febio_spec.Control.step_size = 1/numTimeSteps;
febio_spec.Control.time_stepper.dtmin = dtmin;
febio_spec.Control.time_stepper.dtmax = dtmax;
febio_spec.Control.time_stepper.max_retries = max_retries;
febio_spec.Control.time_stepper.opt_iter = opt_iter;
febio_spec.Control.max_refs = max_refs;
febio_spec.Control.max_ups = max_ups;

%Material section
febio_spec.Material.material{1}.ATTR.type = 'neo-Hookean';
febio_spec.Material.material{1}.ATTR.id = 1;
febio_spec.Material.material{1}.E = youngsMod;
febio_spec.Material.material{1}.v = poisson;

% % % febio_spec.Material.material{2}.ATTR.type='neo-Hookean';
% % % febio_spec.Material.material{2}.ATTR.id=2;
% % % febio_spec.Material.material{2}.E=E_youngs2;
% % % febio_spec.Material.material{2}.v=nu2;

%Geometry section
% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name = 'nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id = (1:size(tibiaMeshV,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL = tibiaMeshV; %The nodel coordinates

% -> Elements
febio_spec.Geometry.Elements{1}.ATTR.type = 'tet4'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat = 1; %material index for this set
febio_spec.Geometry.Elements{1}.ATTR.name = 'tibia'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id = (1:1:size(tibiaMeshE,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL = tibiaMeshE;

% % % febio_spec.Geometry.Elements{2}.ATTR.type='tet4'; %Element type of this set
% % % febio_spec.Geometry.Elements{2}.ATTR.mat=2; %material index for this set
% % % febio_spec.Geometry.Elements{2}.ATTR.name='CancellousBone'; %Name of the element set
% % % febio_spec.Geometry.Elements{2}.elem.ATTR.id=size(E1,1)+(1:1:size(E2,1))'; %Element id's
% % % febio_spec.Geometry.Elements{2}.elem.VAL=E2;

% -> NodeSets
%Top support
febio_spec.Geometry.NodeSet{1}.ATTR.name = 'bcSupportList';
febio_spec.Geometry.NodeSet{1}.node.ATTR.id = bcSupportList(:);
%Force application nodes
febio_spec.Geometry.NodeSet{2}.ATTR.name = 'jointSurfaceNodes';
febio_spec.Geometry.NodeSet{2}.node.ATTR.id = bcForcePrescribeList(:);

%Boundary condition section
% -> Fix boundary conditions
febio_spec.Boundary.fix{1}.ATTR.bc = 'x';
febio_spec.Boundary.fix{1}.ATTR.node_set = febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.fix{2}.ATTR.bc = 'y';
febio_spec.Boundary.fix{2}.ATTR.node_set = febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.fix{3}.ATTR.bc = 'z';
febio_spec.Boundary.fix{3}.ATTR.node_set = febio_spec.Geometry.NodeSet{1}.ATTR.name;
%Force X
febio_spec.MeshData.NodeData{1}.ATTR.name = 'force_X';
febio_spec.MeshData.NodeData{1}.ATTR.node_set = febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.MeshData.NodeData{1}.node.VAL = force_X;
febio_spec.MeshData.NodeData{1}.node.ATTR.lid = (1:1:numel(bcForcePrescribeList))';
%Force Y
febio_spec.MeshData.NodeData{2}.ATTR.name = 'force_Y';
febio_spec.MeshData.NodeData{2}.ATTR.node_set = febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.MeshData.NodeData{2}.node.VAL = force_Y;
febio_spec.MeshData.NodeData{2}.node.ATTR.lid = (1:1:numel(bcForcePrescribeList))';
%Force Z
febio_spec.MeshData.NodeData{3}.ATTR.name='force_Z';
febio_spec.MeshData.NodeData{3}.ATTR.node_set = febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.MeshData.NodeData{3}.node.VAL = force_Z;
febio_spec.MeshData.NodeData{3}.node.ATTR.lid = (1:1:numel(bcForcePrescribeList))';

%Loads section
%Seems like forces are applied consistently in this example
% -> Prescribed nodal forces
febio_spec.Loads.nodal_load{1}.ATTR.bc = 'x';
febio_spec.Loads.nodal_load{1}.ATTR.node_set = febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.Loads.nodal_load{1}.scale.ATTR.lc = 1;
febio_spec.Loads.nodal_load{1}.scale.VAL = 1;
febio_spec.Loads.nodal_load{1}.value.ATTR.node_data = febio_spec.MeshData.NodeData{1}.ATTR.name;

febio_spec.Loads.nodal_load{2}.ATTR.bc = 'y';
febio_spec.Loads.nodal_load{2}.ATTR.node_set = febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.Loads.nodal_load{2}.scale.ATTR.lc = 1;
febio_spec.Loads.nodal_load{2}.scale.VAL = 1;
febio_spec.Loads.nodal_load{2}.value.ATTR.node_data = febio_spec.MeshData.NodeData{2}.ATTR.name;

febio_spec.Loads.nodal_load{3}.ATTR.bc = 'z';
febio_spec.Loads.nodal_load{3}.ATTR.node_set = febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.Loads.nodal_load{3}.scale.ATTR.lc = 1;
febio_spec.Loads.nodal_load{3}.scale.VAL = 1;
febio_spec.Loads.nodal_load{3}.value.ATTR.node_data = febio_spec.MeshData.NodeData{3}.ATTR.name;

%Output section
% -> log file
%Displacement
febio_spec.Output.logfile.ATTR.file = febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file = febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data = 'ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim = ',';
febio_spec.Output.logfile.node_data{1}.VAL = 1:size(tibiaMeshV,1);
%Stress
febio_spec.Output.logfile.element_data{1}.ATTR.file = febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data = 's1;s2;s3';
febio_spec.Output.logfile.element_data{1}.ATTR.delim = ',';
febio_spec.Output.logfile.element_data{1}.VAL = 1:1:size(tibiaMeshE,1); %Rigid body material id
%Strain energy
febio_spec.Output.logfile.element_data{2}.ATTR.file = febioLogFileName_strainEnergy;
febio_spec.Output.logfile.element_data{2}.ATTR.data = 'sed';
febio_spec.Output.logfile.element_data{2}.ATTR.delim = ',';
febio_spec.Output.logfile.element_data{2}.VAL = 1:1:size(tibiaMeshE,1);

%Export to XML file
%Takes a while and likely needs JavaHeap memory increase
%   - Worked pretty well with 2,000Mb allocated to heap (see: https://au.mathworks.com/help/matlab/matlab_external/java-heap-memory-preferences.html)
%Or conversely reduce the number of mesh elements to help here
%This has now been done above to reduce file size (i.e. parts remeshed to
%3.0 edge length). It seems like this is necessary for the simulation to
%not take forever as well.
febioStruct2xml(febio_spec,febioFebFileName);

%% Run the FEBio analysis

%Analysis structure
febioAnalysis.run_filename = febioFebFileName; %The input file name
febioAnalysis.run_logname = febioLogFileName; %The name for the log file
febioAnalysis.disp_on = 1; %Display information on the command window
febioAnalysis.disp_log_on = 1; %Display convergence information in the command window
febioAnalysis.runMode = runMode;%'internal';
febioAnalysis.t_check = 10; %Time for checking log file (dont set too small)
febioAnalysis.maxtpi = 1e99; %Max analysis time
febioAnalysis.maxLogCheckTime = 10; %Max log file checking time

%Run
[runFlag] = runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

%% Import FEBio results
%Note this just comes from the example and could be done in a somewhat more
%efficient manner if desired

if runFlag==1 %i.e. a succesful run
    
    % Importing nodal displacements from log file
    [time_mat, N_disp_mat,~] = importFEBio_logfile(fullfile(pwd,febioLogFileName_disp)); %Nodal displacement
    time_mat = [0; time_mat(:)]; %Time
    
    %Get nodal displacements
    N_disp_mat = N_disp_mat(:,2:end,:);
    sizImport = size(N_disp_mat);
    sizImport(3) = sizImport(3)+1;
    N_disp_mat_n = zeros(sizImport);
    N_disp_mat_n(:,:,2:end) = N_disp_mat;
    N_disp_mat = N_disp_mat_n;
    DN = N_disp_mat(:,:,end);
    DN_magnitude = sqrt(sum(DN(:,3).^2,2));
    V_DEF = N_disp_mat+repmat(tibiaMeshV,[1 1 size(N_disp_mat,3)]);
    
    %Import strain energy density from log file
    [~,E_energy,~] = importFEBio_logfile(fullfile(pwd,febioLogFileName_strainEnergy)); %Element strain energy

    %Remove nodal index column
    E_energy = E_energy(:,2:end,:);

    %Add initial state i.e. zero energy
    sizImport = size(E_energy);
    sizImport(3) = sizImport(3)+1;
    E_energy_mat_n = zeros(sizImport);
    E_energy_mat_n(:,:,2:end) = E_energy;
    E_energy = E_energy_mat_n;
    
    %Convert to appropriate data for visualisation
    [FE_face,C_energy_face] = element2patch(tibiaMeshE,E_energy(:,:,end),'tet4');
    [CV] = faceToVertexMeasure(FE_face,tibiaMeshV,C_energy_face);
    [indBoundary] = tesBoundary(FE_face,tibiaMeshV);
    tibiaMeshFb = FE_face(indBoundary,:);
    
    %Plot the simulated results using anim8
    axLim = [min(min(V_DEF,[],3),[],1); max(max(V_DEF,[],3),[],1)];

    % Create basic view and store graphics handle to initiate animation
    hf = cFigure; %Open figure
    title('Strain energy density')
    gtitle([febioFebFileNamePart,': Press play to animate']);
    hp1 = gpatch(tibiaMeshFb,V_DEF(:,:,end),CV,'k',1); %Add graphics object to animate
    hp1.FaceColor = 'Interp';
    
    %Set colour bar
    axisGeom(gca);
    colormap(gjet(250)); colorbar;
    caxis([0 max(E_energy(:))/25]);
    axis(axLim(:)'); %Set axis limits statically
    camlight headlight;

    % Set up animation features
    animStruct.Time = time_mat; %The time vector
    for qt = 1:1:size(N_disp_mat,3) %Loop over time increments
        
        DN = N_disp_mat(:,:,qt); %Current displacement

        [FE_face,C_energy_face] = element2patch(tibiaMeshE,E_energy(:,:,qt),'tet4');
        [CV] = faceToVertexMeasure(FE_face,tibiaMeshV,C_energy_face);

        %Set entries in animation structure
        animStruct.Handles{qt} = [hp1 hp1]; %Handles of objects to animate
        animStruct.Props{qt} = {'Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt} = {V_DEF(:,:,qt),CV}; %Property values for to set in order to animate
        
    end
    
    %Initiate animation
    anim8(hf,animStruct);

else
    error('FEBio simulation failed...no results to extract...')
end

%% Notes

%%%%% The above simulation seems to work, but it's clear that the set-up,
%%%%% materials etc. aren't great, given that the tibia bends...

