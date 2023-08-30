%-------------------------------------------------------------------------%
% Spacecraft in an asteroid belt
%-------------------------------------------------------------------------%

% Date:    08/2023
% Author:  Ivan Sermanoukian
% Subject:  Bio-inspired Intelligence and learning for Aerospace Applications

%% PREAMBLE

format longE

% Clear workspace, command window and close windows
clear all;
close all;
clc;

% Set LaTeX interpreter
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

% Set the random number generator's seed for reproducibility
rng(42); % You can change the seed value if you'd like

%%
alpha   = 0.6;
gamma   = 0.9;
% eps     = 0.001;
eps = 0;

numTeleports = 15; % Change this to the desired number of teleport tiles

% Possible actions:
UP              = 1;
DOWN            = 2;
RIGHT           = 3;
LEFT            = 4;
FORWARD         = 5;
BACKWARD        = 6;
RIGHT_UP        = 7;
LEFT_DOWN       = 8;
RIGHT_DOWN      = 9;
LEFT_UP         = 10;

actionList      = [UP;DOWN;RIGHT;LEFT;FORWARD;BACKWARD;RIGHT_UP;LEFT_DOWN;RIGHT_DOWN;LEFT_UP];
% actionList      = [UP;DOWN;RIGHT;LEFT;FORWARD;BACKWARD;];

STARTSTATE      = [1,1,1];
GOALSTATE       = [20,20,20];

Nver            = 20;  % Number of states in vertical direction
Nhor            = 20;  % Number of states in horizontal direction
Ndepth          = 20;  % Number of states in horizontal direction

% Generate random positions for teleport tiles
teleportPositions = randi([1, Nver], numTeleports, 1);
teleportCols = randi([1, Nhor], numTeleports, 1);
teleportDepths = randi([1, Ndepth], numTeleports, 1);


% Generate wind perturbations
max_wind = 0;
WIND_VER        = randi([0, 0], 1, Nver);
WIND_HOR        = randi([0, max_wind], 1, Nhor);
WIND_DEPTH      = randi([0, 0], 1, Ndepth);

% Combine the positions into the TELEPORTS matrix
TELEPORTS = [teleportPositions, teleportCols, teleportDepths];

Q               = zeros(Nver,Nhor,Ndepth,size(actionList,1));

n_episodes = 10000;

Policy          = randi(size(actionList,1),Nver,Nhor,Ndepth);
allK            = zeros(n_episodes,1);

% Define blocked layers with specified blocked density for each depth (including depth layer)
blockedLayers = cell(1, Ndepth); % One cell array element per depth layer

middleDensities = rand(Ndepth-2,1)' * 0.6;
blockedDensity = [0, middleDensities, 0]; 

for layer = 1:(Ndepth)
    numBlocked = round(blockedDensity(layer) * Nver * Nhor);
    blockedPositions = randperm(Nver * Nhor, numBlocked);
    [blockedRows, blockedCols] = ind2sub([Nver, Nhor], blockedPositions);
    
    blockedLayer = zeros(Nver, Nhor); % Initialize a blocked layer for this depth
    blockedLayer(sub2ind([Nver, Nhor], blockedRows, blockedCols)) = 1; % Set blocked cells
    blockedLayers{layer} = blockedLayer; % Store blocked layer in the cell array

end

for Episode = 1:n_episodes
    
    % Episode:
    episodeLength       = Nver*Nhor*Ndepth;
    
    % Start at random state:
    state           	= STARTSTATE;
    
    % Action according to current policy:
    actionIdx           = Policy(state(1,1),state(1,2),state(1,3));
    
    action              = actionList(actionIdx);
    
    for k=1:episodeLength
        
        flag_blocked = false;
        if isequal(state, GOALSTATE)
            allK(Episode) = k-1; % reached previous timestep
            break
        end
        
        % Update state based on the action:
        if action == UP
            movementVector = [-1, 0, 0];
        elseif action == DOWN
            movementVector = [1, 0, 0];
        elseif action == RIGHT
            movementVector = [0, 1, 0];
        elseif action == LEFT
            movementVector = [0, -1, 0];
        elseif action == FORWARD
            movementVector = [0, 0, 1];
        elseif action == BACKWARD
            movementVector = [0, 0, -1];
        elseif action == RIGHT_UP
            movementVector = [-1, 1, 0]; 
        elseif action == LEFT_UP
            movementVector = [-1, -1, 0]; 
        elseif action == RIGHT_DOWN
            movementVector = [1, 1, 0];
        elseif action == LEFT_DOWN
            movementVector = [1, -1, 0]; 
        else
            error('Unkown action')
        end
        
        % Calculate the new state:
        newMovedState = state + movementVector;
        newWindState = [1,0,0]*WIND_VER(state(1)) + [0,1,0]*WIND_HOR(state(2)) + [0,0,1]*WIND_DEPTH(state(3));
        newState = newMovedState + newWindState;
        
        if (newState(1) < 1 || newState(1) > Nver) || ...
           (newState(2) < 1 || newState(2) > Nhor) || ...
           (newState(3) < 1 || newState(3) > Ndepth)
            newState = state; % Keep the current state if out of bounds
        end

        % Check if the new state is in a blocked cell in any layer
        if blockedLayers{newState(3)}(newState(1), newState(2))
            newState = state;
            flag_blocked = true;
        else
            % Check if the new state is on a teleport tile
            for current_teleport= 1:1:size(TELEPORTS,1)
                if isequal(newState, TELEPORTS(current_teleport,:)')
                    teleportRow = randi(Nver); % Generate random row
                    teleportCol = randi(Nhor); % Generate random column
                    
                    % Update the agent's position to the teleport destination in the first layer
    %                 newState = [teleportRow, teleportCol, 1];
                    newState = [1, 1, 1];
                end
            end
        end

        % Reward:
        if isequal(newState, GOALSTATE)
            reward = Nver*Nhor*Ndepth;
        elseif flag_blocked == true
            reward=-Nver*Nhor;
        else
            reward=-max(Nhor,Nver);
        end        
        
        % Action according to current policy:
        newActionIdx = Policy(newState(1, 1), newState(1, 2), newState(1, 3));
        newAction = actionList(newActionIdx);
        
        % Update state-action values:
        Q(state(1, 1), state(1, 2), state(1, 3), actionIdx) = Q(state(1, 1), state(1, 2), state(1, 3), actionIdx) + ...
            alpha * (reward + gamma * Q(newState(1, 1), newState(1, 2), newState(1, 3), newActionIdx) - Q(state(1, 1), state(1, 2), state(1, 3), actionIdx));

        % Update Policy:
        if rand < eps
            Policy(state(1, 1), state(1, 2), state(1, 3)) = randi(size(actionList, 1), 1, 1);
        else
            dummy = Q(state(1, 1), state(1, 2), state(1, 3), :);
            [~, Policy(state(1, 1), state(1, 2), state(1, 3))] = max(dummy);
            clear dummy;
        end
        
        % Update state and action:
        state               = newState;
        
        actionIdx           = newActionIdx;
        
        action              = newAction;
        
    end
    if ~isequal(state, GOALSTATE)
        allK(Episode) = episodeLength;
    end
end

%% Policy

figure
scatter(1:1:length(allK),allK,1)

xlabel('Episode', 'FontSize', 18, 'FontWeight', 'bold')
ylabel('Number of steps to goal', 'FontSize', 18, 'FontWeight', 'bold')

%%
% Calculate the aspect ratio for the grid cells
cellWidth = 1;
cellHeight = 1;
aspectRatio = cellWidth / cellHeight;
directionVectors = {[0, -1], [0, 1], [1, 0], [-1, 0], [0, 0], [0, 0], [1, -1], [-1, 1], [1, 1], [-1, -1]}; % Adjust these based on your actionList
% Plot each depth layer separately

newPosition = [1920, 50, 500, 500;
               2420, 50, 500, 500;
               2920, 50, 500, 500;
               1920, 550, 500, 500;
               2420, 550, 500, 500;
                ];

for layer = 1:Ndepth
    fig = figure;
    
    % Create a colormap for the tiles (black for empty, white for blocked)
    tileColormap = [0 0 0; 1 1 1];
    
    imagesc(blockedLayers{layer});
    colormap(flipud(tileColormap));
    caxis([0 1]); % Set colormap limits
    
    % Adjust the aspect ratio of the plot to make each cell a square
    pbaspect([Nhor aspectRatio * Nver 1]);
    
    hold on;
    
    % Plot grid lines to mark tile divisions
    for i = 1:Nver-1
        plot([0.5, Nhor+0.5], [i+0.5, i+0.5], 'k', 'LineWidth', 1); % Horizontal lines
    end
    for j = 1:Nhor-1
        plot([j+0.5, j+0.5], [0.5, Nver+0.5], 'k', 'LineWidth', 1); % Vertical lines
    end

    % Plot "T" for teleport tiles
    for t = 1:size(TELEPORTS, 1)
        teleportRow = TELEPORTS(t, 1);
        teleportCol = TELEPORTS(t, 2);
        if TELEPORTS(t, 3) == layer
            text(teleportCol, teleportRow, 'T', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'green', 'FontSize', 20);
        end
    end
    
    % Plot circles for forward/backward actions, and arrows for other actions
    for i = 1:Nver
        for j = 1:Nhor
            policyActionIdx = Policy(i, j, layer);
            
            % Check if the cell is not blocked
            if blockedLayers{layer}(i, j) == 0 && policyActionIdx ~= 0
                arrowCenter = [j, i]; % Center of the tile
                arrowDirection = directionVectors{policyActionIdx};
                
                if policyActionIdx == FORWARD
                % Add a circle for forward actions
                rectangle('Position', [j-0.4, i-0.4, 0.8, 0.8], 'Curvature', [1 1], 'EdgeColor', 'red', 'LineWidth', 1);
                elseif policyActionIdx == BACKWARD
                    % Add a circle with an "X" for backward actions
                    rectangle('Position', [j-0.4, i-0.4, 0.8, 0.8], 'Curvature', [1 1], 'EdgeColor', 'red', 'LineWidth', 1);
                    text(j, i, 'X', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'red', 'FontSize', 14);
                else
                    % Plot arrows for other actions
                    quiver(arrowCenter(1), arrowCenter(2), arrowDirection(1), arrowDirection(2), 'Color', 'red', 'LineWidth', 1, 'AutoScale', 'off','MaxHeadSize', 0.8);
                end
            end
            
            % Mark initial and final tiles with labels
            if [i, j, layer] == STARTSTATE
                text(j, i, 'S', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'blue', 'FontSize', 20);
            elseif [i, j, layer] == GOALSTATE
                text(j, i, 'E', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'blue', 'FontSize', 20);
            end
        end
    end
    
    hold off;
    
    xlabel('Horizontal');
    ylabel('Vertical');
    title(['Policy Arrows - Depth ' num2str(layer)]);
    
    % Add colorbar to indicate tile types
    colorbar('Ticks', [0.25, 0.75], 'TickLabels', {'Empty', 'Blocked'});
end




