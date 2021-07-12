

cd('..\ShapeModel\shape_model');

%Load mean STL
[tibiaSTLstruct] = import_STL('tibia_mean.stl');
meanF = tibiaSTLstruct.solidFaces{1}; %Faces
meanV = tibiaSTLstruct.solidVertices{1}; %Vertices
[meanF,meanV] = mergeVertices(meanF,meanV);

% % % %Visualise
% % % cFigure;
% % % gpatch(meanF,meanV);
% % % axisGeom;

%Read in created points from gias2 reconstruction
meanPoints = readmatrix('pointCloud_mean.csv');
for cc = 1:10
    
    %Load the reconstructed component points
    p3_Points = readmatrix(['pointCloud_pc',num2str(cc),'_p3.csv']);
    m3_Points = readmatrix(['pointCloud_pc',num2str(cc),'_m3.csv']);

    %%%%%%%% TODO: better looping...

    %The reconstructed points don't align with the STL points after merging,
    %but theoretically we should be able to find the matching points based on
    %the minimum distance. We calculate the distance between each individual
    %point in the mean imported point cloud and take the minimum to be the
    %matching point. This gives us the index of the row to match the data to
    %across the mean and reconstructed points so that they then match the STL.
    for pp = 1:length(meanPoints)

        %Calculate current point distances to STL nodes
        distCurrPoint = distancePoints3d(meanV(pp,:),meanPoints);

        %Get the index of the min distance points
        minInd = find(distCurrPoint == min(distCurrPoint));

        %Place the points from the reconstructed point clouds at the
        %appropriate row to match the STL
        %Note that we don't really need to reallign the mean points here as
        %they already exist from the STL. We can check it here for sanity
        %though if we want...
    % % %     meanCheck(pp,:) = meanPoints(minInd,:);
        pcPlusPoints(pp,:) = p3_Points(minInd,:);
        pcMinusPoints(pp,:) = m3_Points(minInd,:);

    end
    clear pp

% % %     %Visualise
% % %     cFigure;
% % %     gpatch(meanF,pcPlusPoints);
% % %     axisGeom;

    %Calculate distances between points for heat mapping
    pcPlusDist = distancePoints3d(meanV,pcPlusPoints);
    pcMinusDist = distancePoints3d(meanV,pcMinusPoints);

    %Convert distances to face measures
    pcPlusC = vertexToFaceMeasure(meanF,pcPlusDist);
    pcMinusC = vertexToFaceMeasure(meanF,pcMinusDist);

    %Visualise surfaces and change in points

    %Create figure
    cFigure;

    %Minus PC
    subplot(1,4,1);
    title(['PC',num2str(cc),' -3 SD'],'fontsize',10);
    gpatch(meanF,pcMinusPoints,'rw','none')
    axisGeom; view(2); axis off
    camlight headlight

    %Mean
    subplot(1,4,2);
    title('Mean','fontsize',10);
    gpatch(meanF,meanV,'bw','none')
    axisGeom; view(2); axis off
    camlight headlight

    %Plus PC
    subplot(1,4,3);
    title(['PC',num2str(cc),' +3 SD'],'fontsize',10);
    gpatch(meanF,pcPlusPoints,'gw','none')
    axisGeom; view(2); axis off
    camlight headlight

    %Heatmap
    subplot(1,4,4);
    title(['Change Map PC',num2str(cc)]);
    gpatch(meanF,pcPlusPoints,pcPlusC,'none');
    colormap viridis; %colorbar;
    axisGeom; view(2); axis off;
    camlight headlight

    %Reset all axes so that they are the same
    %Get the limits from each axes
    for aa = 1:4
        subplot(1,4,aa);
        ax = gca;
        xLim(aa,:) = ax.XLim;
        yLim(aa,:) = ax.YLim;
        zLim(aa,:) = ax.ZLim;    
    end
    clear aa
    %Get the min and max from each
    setX(1,1) = min(xLim(:,1)); setX(1,2) = max(xLim(:,2));
    setY(1,1) = min(yLim(:,1)); setY(1,2) = max(yLim(:,2));
    setZ(1,1) = min(zLim(:,1)); setZ(1,2) = max(zLim(:,2));
    %Set the axes values
    for aa = 1:4
        subplot(1,4,aa);
        ax = gca;
        ax.XLim = setX;
        ax.YLim = setY;
        ax.ZLim = setZ;
    end
    clear aa
    
end
clear cc



