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

% Load data
load("data_10x10x5_100000_500_25.mat")

% Set limits
global_allk = global_allk(:,64:1260);
global_time = global_time(64:1260);
variables = variables(:,64:1260);

%%
min_vec = zeros(size(global_allk,2),1);
episodes_vec = zeros(size(global_allk,2),1);
alpha = length(variables(1,:));
gamma = length(variables(1,:));
eps = length(variables(1,:));
for i=1:1:size(global_allk,2)

    alpha(i)   = variables(1,i);
    gamma(i)   = variables(2,i);
    eps(i)     = variables(3,i);

    values = global_allk(global_allk(:,i) > 0,i);
    min_vec(i) = min(values(values > 0));
    mode_vec(i) = mode(values(values > 0));
    episodes_vec(i) = length(values);
end

%%
fig1 = figure(1)
alpha_lin = linspace(min(alpha), max(alpha), 75)
gamma_lin = linspace(min(gamma), max(gamma), 75)
[X,Y] = meshgrid(alpha_lin, gamma_lin);
Z = griddata(alpha,gamma,min_vec,X,Y,'v4');
mesh(X,Y,Z)
axis tight; hold on
colorbar;
caxis([12 17]);
xlim([0.05 1])
ylim([0.1 1])
zlim([12 17])
xlabel("$\alpha$")
ylabel("$\gamma$")
zlabel("Minimum number of steps")
fontsize(fig1, 16, "points")

%%
fig2 = figure(2)
Z = griddata(alpha,gamma,full(mode_vec),X,Y,'v4');
mesh(X,Y,Z)
axis tight; hold on
xlim([0.05 0.95])
ylim([0.05 1])
zlim([12 20])
colorbar;
caxis([12 20]);
xlabel("$\alpha$")
ylabel("$\gamma$")
zlabel("Modal number of steps")
fontsize(fig2, 16, "points")

%%
fig4 = figure(4)
alpha_lin = linspace(min(alpha), max(alpha), 75);
gamma_lin = linspace(min(gamma), max(gamma), 75);
[X,Y] = meshgrid(alpha_lin, gamma_lin);
Z = griddata(alpha,gamma,global_time,X,Y,'v4');
mesh(X,Y,Z)
axis tight; hold on
colorbar;
caxis([15 250]);
% 
xlabel("$\alpha$")
ylabel("$\gamma$")
zlabel("Global time [s]")
fontsize(fig4, 16, "points")

