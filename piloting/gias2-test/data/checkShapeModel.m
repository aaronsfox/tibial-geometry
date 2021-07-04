%Navigate to shape model file
cd('shape_model')

%Load the mean stl
[meanSTLstruct] = import_STL('femur_mean.stl');
meanF = meanSTLstruct.solidFaces{1}; %Faces
meanV = meanSTLstruct.solidVertices{1}; %Vertices
[meanF,meanV] = mergeVertices(meanF,meanV);

%Visualise
h = cFigure; hold on
gpatch(meanF,meanV,'gw','k');
axisGeom;

%Load +/- SD for PC1
[plusSTLstruct] = import_STL('femur_recon_pc1p2.stl');
plusF = plusSTLstruct.solidFaces{1}; %Faces
plusV = plusSTLstruct.solidVertices{1}; %Vertices
[plusF,plusV] = mergeVertices(plusF,plusV);

[minusSTLstruct] = import_STL('femur_recon_pc1m2.stl');
minusF = minusSTLstruct.solidFaces{1}; %Faces
minusV = minusSTLstruct.solidVertices{1}; %Vertices
[minusF,minusV] = mergeVertices(minusF,minusV);

%Add to visualisation
gpatch(plusF,plusV,'rw','k');
gpatch(minusF,minusV,'bw','k');

%% Compare min distance algorithm between pc

%plus PC1 to mean

%Compare distances
% % % [D2] = triSurfSetDist(plusF,plusV,meanF,meanV,'dist');
[D2] = triSurfSetDist(minusF,minusV,meanF,meanV,'dist');

%Visualise
% % % [CF] = vertexToFaceMeasure(plusF,D2);
[CF] = vertexToFaceMeasure(minusF,D2);

hf2=cFigure;
% title('Closest point distance metric on plus PC1');
hold on;
% % % patch('faces',plusF,'vertices',plusV,'FaceColor','flat','CData',CF);
patch('faces',minusF,'vertices',minusV,'FaceColor','flat','CData',CF);
% % % patch('faces',F1,'vertices',V1,'FaceColor',0.5.*ones(1,3),'FaceAlpha',faceAlpha1,'EdgeColor','None');
gpatch(meanF,meanV,'k','none',0.5)

colormap jet; colorbar;
axis equal; view(3); axis tight; axis off;
set(gca,'FontSize',font_size);
camlight headlight;
drawnow;


