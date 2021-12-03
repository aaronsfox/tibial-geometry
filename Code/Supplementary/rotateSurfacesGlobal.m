function [tibiaF, tibiaV, fibulaF, fibulaV, trabF, trabV] = rotateSurfacesGlobal(landmarks, noFib, tibiaF, tibiaV, fibulaF, fibulaV, trabF, trabV)

    %% Function that simplifies/cleans up the main code for rotating surfaces
    %  based on the extracted landmarks
    
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
    for pp = 1:length(trabV)
        trabV(pp,:) = transformPoint3d(trabV(pp,:),globalTransform);    
    end
    clear pp
    if ~noFib
        for pp = 1:length(fibulaV)
            fibulaV(pp,:) = transformPoint3d(fibulaV(pp,:),globalTransform);    
        end
        clear pp
    end

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
    for pp = 1:length(trabV)
        trabV(pp,:) = transformPoint3d(trabV(pp,:),rotZ);    
    end
    clear pp
    if ~noFib
        for pp = 1:length(fibulaV)
            fibulaV(pp,:) = transformPoint3d(fibulaV(pp,:),rotZ);    
        end
        clear pp
    end

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
    for pp = 1:length(trabV)
        trabV(pp,:) = transformPoint3d(trabV(pp,:),transMatrix);    
    end
    clear pp
    if ~noFib
        for pp = 1:length(fibulaV)
            fibulaV(pp,:) = transformPoint3d(fibulaV(pp,:),transMatrix);    
        end
        clear pp
    end

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
    for pp = 1:length(trabV)
        trabV(pp,:) = transformPoint3d(trabV(pp,:),createRotationOz(deg2rad(-90)));    
    end
    clear pp
    if ~noFib
        for pp = 1:length(fibulaV)
            fibulaV(pp,:) = transformPoint3d(fibulaV(pp,:),createRotationOz(deg2rad(-90)));    
        end
        clear pp
    end

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
    for pp = 1:length(trabV)
        trabV(pp,:) = transformPoint3d(trabV(pp,:),createRotationOy(deg2rad(-90)));    
    end
    clear pp
    if ~noFib
        for pp = 1:length(fibulaV)
            fibulaV(pp,:) = transformPoint3d(fibulaV(pp,:),createRotationOy(deg2rad(-90)));    
        end
        clear pp
    end

    %Transform landmarks
    currLandmarks = fieldnames(landmarks);
    for ff = 1:length(currLandmarks)
        landmarks.(currLandmarks{ff}) = transformPoint3d(landmarks.(currLandmarks{ff}),createRotationOy(deg2rad(-90)));
    end
    clear ff

end
