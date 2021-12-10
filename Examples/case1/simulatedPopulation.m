%% Case 1 - Simulated populations 
%% Step 1 -  Load data create seed 

load('tibiaShapeModel.mat');

%% Step two - determine inputs 

%Number of simulated structures 
nS = 12;

%SD bounds 

SDbounds = 3;

%PC variables to be used 
% Tibia only 
PCselect = [true, true, true, true, true];

% Tibia and Fibula
%%%PCselect = [true, false, false, false, false, false, false];


%% Step 3 - Generate random scores for each PC 
rng(13)
%generate random numbers between SDbounds
for PCind = 1:length(PCselect)
    
    if PCselect(1, PCind)
       r = (SDbounds*rand(nS,PCind));  
    else 
       r(nS,PCind) = zeros;
    end 
    
end
%Convert Numbers to raw score using SD
tibiaSimulatedModels.rawScore = r .* std(tibiaShapeModel.score(:,1:tibiaShapeModel.retainPCs));

%I think this is the reconstruction step? Maybe? 
%may need to rename
tibiaSimulatedModels.reconstructedScore = tibiaSimulatedModels.rawScore * tibiaShapeModel.loadings(:,1:tibiaShapeModel.retainPCs)';

%Add mean to each simulated tibia 
tibiaSimulatedModels.reconstructedMean = tibiaSimulatedModels.reconstructedScore (1:nS,:) + tibiaShapeModel.mean;


%% Visualisation 

%Set array to rotate views of surfaces later in reconstruction process
surfaceRot = [-90, 0, 90, 180];

    
%loop through each simulated tibia to create subplot 
    
     cFigure; hold on;
     
for nSind = 1:nS
    
    %create subplot
        subplot(2,6,nSind); hold on;
        
    %Add surface to subplot 
        hp = gpatch(tibiaShapeModel.F,reshape(tibiaSimulatedModels.reconstructedMean(nSind,:),...
            [3, length(tibiaSimulatedModels.reconstructedMean)/3])','#e3dac9','none', 0.5);
        
        %Set axis view
        axis equal; axis tight;
        view(0,90); rotate(hp,[0 1 0], surfaceRot(1));
        
        %Set axis parameters
        camlight headlight; axis on
        
        %Add title
        title([{nSind}], 'FontSize', 12);
    
end

export_fig(['testfig-popsize12-allPC-3SD.png']);
	close