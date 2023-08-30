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
alpha_vec   = linspace(0,1,11);
gamma_vec   = linspace(0,1,11);
eps_vec     = linspace(0,0.01,2);

n_episodes = 10000;
flag_max = 500;
variables = zeros(3,length(alpha_vec)*length(gamma_vec)*length(eps_vec));
global_allk = sparse(n_episodes,length(alpha_vec)*length(gamma_vec)*length(eps_vec));
global_time = zeros(length(alpha_vec)*length(gamma_vec)*length(eps_vec),1);
%%
numTeleports = 25; % Change this to the desired number of teleport tiles

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
GOALSTATE       = [5,5,5];

Nver            = 5;  % Number of states in vertical direction
Nhor            = 5;  % Number of states in horizontal direction
Ndepth          = 5;  % Number of states in horizontal direction

% Generate random positions for teleport tiles
teleportPositions = randi([1, Nver], numTeleports, 1);
teleportCols = randi([1, Nhor], numTeleports, 1);
teleportDepths = randi([1, Ndepth], numTeleports, 1);


% Generate wind perturbations
max_wind = 0;
WIND_VER        = randi([0, max_wind], 1, Nver);
WIND_HOR        = randi([0, max_wind], 1, Nhor);
WIND_DEPTH      = randi([0, max_wind], 1, Ndepth);

% Combine the positions into the TELEPORTS matrix
TELEPORTS = [teleportPositions, teleportCols, teleportDepths];

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

%%
counter = 0;
for alpha_index = 1:1:length(alpha_vec)
    alpha = alpha_vec(alpha_index);
    for gamma_index = 1:1:length(gamma_vec)
        gamma = gamma_vec(gamma_index);
        for eps_index = 1:1:length(eps_vec)
            tic;
            eps = eps_vec(eps_index);
            counter = counter + 1;
            variables(1,counter) = alpha;
            variables(2,counter) = gamma;
            variables(3,counter) = eps;

            %
            Q               = zeros(Nver,Nhor,Ndepth,size(actionList,1));
            
            Policy          = randi(size(actionList,1),Nver,Nhor,Ndepth);
            allK            = zeros(n_episodes,1);
            
            Episode = 1;
            flag = 1;
            while (Episode <= n_episodes) && (flag <= flag_max)
                
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
                if Episode > 1
                    if (allK(Episode) == allK(Episode - 1))
                        flag = flag + 1;
                    else
                        flag = 1;
                    end
                end
                Episode = Episode + 1;
            end
            global_allk(:,counter) = allK;
            global_time(counter) = toc;

        end
    end
end

save("data.mat")




