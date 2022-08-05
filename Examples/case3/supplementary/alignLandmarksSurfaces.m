function [imageLandmarks, tibiaV, fibulaV] = alignLandmarksSurfaces(imageLandmarks, tibiaV, fibulaV)

    %% Function that simplifies/cleans up the main code for aligning landmarks 
    %  and surfaces to the global coordinate system
    
    %% Image data - surfaces and landmarks
    
    %Create the tibial planes
    %Frontal
    planes.frontal = createPlane(imageLandmarks.IM,...
        imageLandmarks.LC,imageLandmarks.MC);
    %Torsional
    planes.torsional = createPlane(imageLandmarks.IC,...
        imageLandmarks.MM,imageLandmarks.LM);

    %Create transform to get tibia aligned to the global plane
    globalTransform = createBasisTransform3d('global',planes.torsional);

    %Transform surfaces
    tibiaV = transformPoint3d(tibiaV,globalTransform);
    fibulaV = transformPoint3d(fibulaV,globalTransform);
    
    %Transform landmarks
    currLandmarks = fieldnames(imageLandmarks);
    for ff = 1:length(currLandmarks)
        imageLandmarks.(currLandmarks{ff}) = transformPoint3d(imageLandmarks.(currLandmarks{ff}),globalTransform);
    end
    clear ff

    %Do the secondary rotation around the X-axis to make the tibia vertical
    %along the Y-axis

    %Identify distance between IC and IM, along XY plane
    pt1 = [imageLandmarks.IC(1),imageLandmarks.IC(2)];
    pt2 = [imageLandmarks.IM(1),imageLandmarks.IM(2)];
    dIC_IM = sqrt((pt2(2) - pt1(2))^2 + (pt2(1) - pt1(1))^2);

    %Identify distance of IC from IM along the x-axis
    pt3 = [imageLandmarks.IC(1),imageLandmarks.IM(2)];
    dIC_IMx = sqrt((pt3(2) - pt1(2))^2 + (pt3(1) - pt1(1))^2);

    %Calculate angle to rotate about x-axis
    rotAng = asin(dIC_IMx/dIC_IM);

    %Create rotation matrix around Z-axis by specified angle (in radians)
    rotZ = createRotationOz(rotAng*-1); %-ve for anti-clockwise

    %Transform surfaces
    tibiaV = transformPoint3d(tibiaV,rotZ);    
    fibulaV = transformPoint3d(fibulaV,rotZ);    

    %Transform landmarks
    for ff = 1:length(currLandmarks)
        imageLandmarks.(currLandmarks{ff}) = transformPoint3d(imageLandmarks.(currLandmarks{ff}),rotZ);
    end
    clear ff

    %Create the transform to make the IM landmark the global origin
    transMatrix = createTranslation3d([0,0,0] - imageLandmarks.IM);

    %Transform surfaces
    tibiaV = transformPoint3d(tibiaV,transMatrix);    
    fibulaV = transformPoint3d(fibulaV,transMatrix);    
    
    %Transform landmarks
    for ff = 1:length(currLandmarks)
        imageLandmarks.(currLandmarks{ff}) = transformPoint3d(imageLandmarks.(currLandmarks{ff}),transMatrix);
    end
    clear ff

    %The current bodies are aligned so that X is vertical and Y is lateral
    %This needs to be shifted to align with the ISB recommendations so that Y
    %is vertical and Z is lateral. This can easily be done by a few rotations
    %about specific axes

    %First, rotate about the z-axis by -90 degrees

    %Transform surfaces
    tibiaV = transformPoint3d(tibiaV,createRotationOz(deg2rad(-90)));    
    fibulaV = transformPoint3d(fibulaV,createRotationOz(deg2rad(-90)));    
    
    %Transform landmarks
    for ff = 1:length(currLandmarks)
        imageLandmarks.(currLandmarks{ff}) = transformPoint3d(imageLandmarks.(currLandmarks{ff}),createRotationOz(deg2rad(-90)));
    end
    clear ff

    %Second, rotate about the y-axis by -90 degrees

    %Transform surfaces
    tibiaV = transformPoint3d(tibiaV,createRotationOy(deg2rad(-90)));    
    fibulaV = transformPoint3d(fibulaV,createRotationOy(deg2rad(-90)));

    %Transform landmarks
    for ff = 1:length(currLandmarks)
        imageLandmarks.(currLandmarks{ff}) = transformPoint3d(imageLandmarks.(currLandmarks{ff}),createRotationOy(deg2rad(-90)));
    end
    clear ff
    
% % %     %Visualise to check
% % %     cFigure; hold on;
% % %     %Tibia
% % %     plotV(tibiaV, 'g.')
% % %     %Fibula
% % %     plotV(fibulaV, 'b.')
% % %     %Image landmarks
% % %     currLandmarks = fieldnames(imageLandmarks);
% % %     for ff = 1:length(currLandmarks)
% % %         plotV(imageLandmarks.(currLandmarks{ff}), 'r.', 'MarkerSize', 25)
% % %     end
% % %     %Axis parameters
% % %     axisGeom;

end
