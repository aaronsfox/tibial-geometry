function [reconstructed] = simulatedPopulation(nSamples, shapeModel, selectPCs, sdBounds, opts)

    %% This function allows users to generate a selected number of surfaces
    %  from a pre-loaded shape model. Note that the below example focuses on
    %  the tibia shape model only, however this process could be applied to
    %  any of the shape models included in this repository. See the
    %  README.MD in this folder for more descriptive details on this
    %  process.
    %
    %  Authors:
    %     Meghan Keast
    %     Centre for Sport Research
    %     Deakin University
    %     mfkeast@deakin.edu.au
    %
    %     Aaron Fox
    %     Centre for Sport Research
    %     Deakin University
    %     aaron.f@deakin.edu.au
    %
    %   Function inputs:
    %       nSamples        the desired number of sample surfaces to generate
    %       shapeModel      a structure with the pre-loaded shape mode
    %       selectPCs       a boolean array matching up to the PCs to
    %                       reconstruct the surface from. This can be a
    %                       maximum length of the number of total PCs
    %                       available in the model
    %       sdBounds        a single value for the standard deviation bounds
    %                       to manipulate the PCs within that will be
    %                       applied in a +/- value (i.e. sdBounds = 3 will
    %                       result in bounds of +/- 3SD) (default = 2)
    %       opts            options structure containing optional inputs of:
    %                       > seed: value to set random seed for consistency (default = none)
    %                       > exportSTL: true/false whether to export surfaces as STL to /stl directory (default = false)
    %                       > exportIMG: true/false whether to export an image grid of the resulting surfaces (default = false)
    %                       > heatMapIMG: true/false whether to add heatmap colouring relative to mean for exported image (default = false)
    %                       > labelsOnIMG: true/false whether to label surfaces with numbers on image (default = true) 
    %
    %   Function outputs:
    %       reconstructed   an n x 3 x nSamples array (where n = the number
    %                       of nodes in the shape model) of the
    %                       reconstructed points from the process
    %
    %% EXAMPLE
    %   
    %       %Load the tibia and trabecular shape models required
    %       load('..\..\ShapeModels\tibia\tibiaShapeModel.mat');
    %
    %       %Set the desired number of samples
    %       nSamples = 30;
    %
    %       %Set the PCs to use in the reconstruction
    %       %This example uses the first 5, ignores the 6th, and uses the 7th
    %       selectPCs = [true, true, true, true, true, false, true];
    %
    %       %Set the standard deviation bounds
    %       sdBounds = 1.5;
    %
    %       %Set the options
    %       opts.seed = 12345;
    %       opts.exportSTL = true;
    %       opts.exportIMG = true;
    %       opts.heatMapIMG = true;
    %       opts.labelsOnIMG = false;
    %
    %       %Run function
    %       simulatedPopulation(nSamples, tibiaShapeModel, selectPCs, sdBounds, opts)

    %% Check function inputs

    %Check that there are at least 3 inputs
    if nargin < 3
        error('At least the number of samples, shape model structure and selected PCs are needed for function.');
    end

    %Check that shape model is a struct
    if ~isstruct(shapeModel)
        error('Shape model input must be a structure containing the shape model data')
    else
        %Check that the shape model structure has the necessary fields
        if ~isfield(shapeModel,'loadings')
            error('Shape model structure must have a field representing the coefficient loadings.')
        end
        if ~isfield(shapeModel,'score')
            error('Shape model structure must have a field representing the shape model sample scores.')
        end
        if ~isfield(shapeModel,'mean')
            error('Shape model structure must have a field representing the shape model nodes mean.')
        end
        if opts.exportSTL || opts.exportIMG
            if ~isfield(shapeModel,'F')
                error('Shape model structure must have representative faces to generate STL and/or image.')
            end
        end
        if opts.exportIMG && opts.heatMapIMG
            if ~isfield(shapeModel,'meanPoints')
                error('Shape model structure must have n x 3 mean points to generate heatmap image.')
            end
        end
    end

    %Check that selectPCs isn't longer than total PCs
    if length(selectPCs) > size(shapeModel.loadings,2)
        error('Number of PCs to generate surfaces from is greater than the number of PCs in the shape model.')
    end

    %Set defaults for sd bounds if necessary
    if nargin < 4
        %Set SD bounds to default
        sdBounds = 2;
    else
        %Check if sdBounds is a single numerical value
        if ~(isnumeric(sdBounds) && length(sdBounds) == 1)
            error('SD bounds must be a single numerical value.')
        end
    end

    %Check for options structure
    if nargin < 5
        %Set default for all options
        opts.exportSTL = false;
        opts.exportIMG = false;
        opts.heatMapIMG = false;
        opts.labelsOnIMG = true;
    else
        %Check for fields within options and set defaults otherwise
        if ~isfield(opts,'exportSTL')
            opts.exportSTL = false;
        end
        if ~isfield(opts,'exportIMG')
            opts.exportIMG = false;
        end
        if ~isfield(opts,'heatMapIMG')
            opts.heatMapIMG = false;
        end
        if ~isfield(opts,'labelsOnIMG')
            opts.labelsOnIMG = true;
        end
    end
    
    %% Set-up
    
    %Turn off warnings
    warning off

    %% Generate random scores for each selected PC 
    
    %Set the seed if desired
    if isfield(opts,'seed')   
        rng(opts.seed)
    end
    
    %Generate nSamples random values from a uniform distribution
    %between the desired SD bounds
    for PCind = 1:length(selectPCs)
        
        %Check if current PC is being used in sampling
        if selectPCs(PCind)
            
            %Get the random numbers between the specified SD bounds
            r(:,PCind) = -sdBounds + (sdBounds - -sdBounds) * rand(nSamples,1);
            
        else
            
            %Set to zeros
            r(:,PCind) = zeros(nSamples,1);

        end
    
    end
    
    %Convert SD values to raw score using SD from the shape model
    simulatedScores = r .* std(shapeModel.score(:,1:length(selectPCs)));
    
    %Reconstruct each case using the simulated PC scores
    simulatedNodes = simulatedScores(:,1:length(selectPCs)) * ...
        shapeModel.loadings(:,1:length(selectPCs))';
    
    %Reshape and add the mean
    for sampleInd = 1:nSamples
        reconstructed(:,:,sampleInd) = reshape(simulatedNodes(sampleInd,:) + ...
            shapeModel.mean, [3, length(shapeModel.mean)/3])';
    end

    %% Export STLs
    
    if opts.exportSTL
        
        %Create a directory to store the STL files in
        mkdir('stl');
        
        %Setup the waitbar
        wbar = waitbar(0,'Exporting simulated surfaces....');
        
        %Loop through the samples
        for sampleInd = 1:nSamples
            
            %Create stl structure
            stlStruct.solidNames = {['simulatedSurface_',num2str(sampleInd)]}; %names of parts
            stlStruct.solidVertices = {reconstructed(:,:,sampleInd)}; %Vertices
            stlStruct.solidFaces = {shapeModel.F}; %Faces
            stlStruct.solidNormals={[]};
            
            %Export the STL
            export_STL_txt(['stl\simulatedSurface_',num2str(sampleInd),'.stl'], stlStruct);

            %Update waitbar
            waitbar(sampleInd/nSamples, wbar);
            
        end
        
        %Update and close waitbar
        waitbar(1, wbar, 'Exported all surfaces. Finishing...');
        pause(1); close(wbar);
       
    end
    
    %% Export image
    
    if opts.exportIMG
        
        %Identify the maximum and minimum axes values based on the minimum
        %and maximum values along the visible X and Y axes
        xMin = min(min(reconstructed(:,1,:)));
        xMax = max(max(reconstructed(:,1,:)));
        yMin = min(min(reconstructed(:,2,:)));
        yMax = max(max(reconstructed(:,2,:)));
        
        %Determine a relatively even square grid to subplot onto based on
        %the number of simulated samples
        subplotSq = sqrt(nSamples);
        
        %Create figure
        cFigure; hold on;

        %Loop through each simulated tibia to create subplot 
        for sampleInd = 1:nSamples
            
            %create subplot
            subplot(ceil(subplotSq), round(subplotSq), sampleInd); hold on;
            
            %Add surface to subplot
            if opts.heatMapIMG
                
                %Calculate point error distance
                errorDist = distancePoints3d(shapeModel.meanPoints, reconstructed(:,:,sampleInd));
                
                %Convert distance error to colour scales for visualisation
                errorDistColF = vertexToFaceMeasure(shapeModel.F, errorDist);
                errorDistColV = faceToVertexMeasure(shapeModel.F, shapeModel.meanPoints, errorDistColF);
                
                %Add surface with heatmap colouring
                hp = gpatch(shapeModel.F, reconstructed(:,:,sampleInd), ...
                    errorDistColV, 'none', 1);

                %Interpolate colouring for smoothness
                hp.FaceColor = 'Interp'; colormap viridis
                                
            else
                
                %Add surface with bone colouring
                hp = gpatch(shapeModel.F, reconstructed(:,:,sampleInd), ...
                    '#e3dac9', 'none', 1);
                
            end
            
            %Set axis view
            axis equal; axis tight;
            view(0,90); 
            rotate(hp,[0 1 0], -90)
            camlight headlight
            
            %Set axis limits
            ax = gca;
            ax.XLim = [xMin xMax];
            ax.YLim = [yMin yMax];
            
            %Turn axis off
            axis off;

            %Add title if required
            if opts.labelsOnIMG
                title(num2str(sampleInd), 'FontSize', 250 / nSamples);
            end
                
        end
        
        %Export figure
        mkdir('img');
        export_fig('img\simulatedPopulation.png','-m5');
        close all

    end
    
%% ----- End of simulatedPopulation.m ----- %%
end