%------------------------------------------------------------------
%                  W5. The final model
%------------------------------------------------------------------

% Parameters:
param.D = 40; % Diffusivity (m2/d) Huismann (43 m2/d)
param.u = 1; % Settling velocity (m/day) Huisman (0.96 m/d)

param.I0 = 350*(60*60*24)/1000; % Light intensity at surface (mmol photons/m^2/day) Huismann (350 µmol phot/m2/s)
% dividing with 1000 to convert from µmol --> mmol
param.k = 15e-12; % spec. light attenuation of phytoplankton (m^2/(mmol N)), HUismann (15e-12 m2/cell)
param.Kbg = 0.2; % light attenuation by sea water/background turbidity (1/m), Huismann (0.2 1/m)

param.H_I = 30*(60*60*24)/1000; % half-saturation constant of light-limited growth (mmol photons/m^2/day), Huismann (30 µmol phot/m2/s)
param.gmax = 0.96; % maximum growth rate (1/day) Huismann (0.96 1/d), Huismann (0.96 1/d)
param.m = 0.24; % loss rate (1/day), Huismann (0.01 1/h)

param.alpha = 1e-9; % nutrient content of phytoplankton cell (mmol N/cell) 
param.eps = 0.5; % phytoplankton recycling coefficient/remineralization rate (dimensionless)
param.H_N = 0.0425; % half-saturation constant of nutrient-limited growth (mmol N/m^2/day)
param.nz = 50; % nutrient concentration at bottom (mmol N/m^3)


param.depth = 100; % water column depth (m)
param.n = 50; %no of grid cells
param.dz = param.depth/param.n; % width of grid cells/grid spacing (m)
param.z = param.dz/2:param.dz:(param.depth-param.dz/2);


% Initial conditions
P0 = 50e7*exp(-(param.z-100).^2/6);
N0 = ones(1,param.n);

% Run model without seasonal forcing
options = odeset('nonnegative', param.n*1:3);
tic
[t, y] = ode45(@PmodelDerivNP, [0,800], [P0,N0], options, param);
toc
P = y(:,1:param.n); 
N = y(:,param.n+1:end);
I = (CalcLight(P(end,:),param))';

%%  Plotting results
% Surface plot
clf
subplot(2,1,1)
surface(t, -param.z,P')
shading flat
xlabel('time (days)')
ylabel('depth (m)')
title('Phytoplankton (mmol N m^-^3)')
colorbar

subplot(2,1,2)
surface(t, -param.z,N')
shading flat
xlabel('time (days)')
ylabel('depth (m)')
title('Nutrients (mmol N m^-^3)')
colorbar

%% Grid sensitivity plot
% clf
% ngrid = [20 50 100 200 300];
% for i = 1:length(ngrid)
%     param.n = ngrid(i);
%     [t, y] = ode45(@PmodelDerivNP, [0,800], [P0,N0], options, param);
%         P = y(:,1:param.n);
%         N = y(:,param.n+1:end);
%
%     plot(P(end,:), -param.z,'o-')
%     drawnow
%     hold on
%     %plot(N(end,:),-param.z,'b','x-')
% end
% xlabel('Concentration (mmol N m^-^3)')
% ylabel('Depth (m)')

%% 3 profile plots

% calculate limiting factor:
I_lim = I./(param.H_I + I);
N_lim = N(end,:)./(param.H_N + N(end,:));
I_lim2 = I_lim

clf
subplot(1,3,1)
plot(P(end,:),-param.z, 'g', 'linewidth',1)
hold on
plot(N(end,:)*1.5e8,-param.z, 'b', 'linewidth',1)
xlabel('mmol N m^-^3')
ylabel('Depth (m)')

subplot(1,3,2)
plot(I,-param.z,'black','linewidth',1)
xlabel('mmol photons m^2 d^-^1')

subplot(1,3,3)
plot(I_lim, -param.z,'color','#EDB120','linewidth',1)
hold on 
plot(N_lim, -param.z,'color','#0072BD','linewidth',1)
hold on
yline(-4,'--','depth = 4 m','linewidth',1)
xlabel('Limiting factor')

%% Sensitivity analysis gmax=µmax
% For-loop
clf
gmax = 0.5:0.1:1;
for i = 1:length(gmax)
    param.gmax = gmax(i);
    [t, y] = ode45(@PmodelDerivNP, [0,800], [P0,N0], options, param);
        P = y(:,1:param.n);
        N = y(:,param.n+1:end);

        plot(P(end,:),-param.z, 'g', 'linewidth',i/2);
        hold on
        plot(N(end,:)*1.5e8,-param.z, 'b', 'linewidth',i/2);
end
title('Sensitivity analysis of µmax')
xlabel('mmol N/m^-^3')
ylabel('Depth (m)')
legend('phytoplankton','nutrients','','','','','','','','','','')

%% Seasonality
% Run model WITH seasonal forcing
tic
[t, y] = ode45(@PmodelDerivNPSeason, [0,2000], [P0,N0], options, param);
toc
P = y(:,1:param.n); 
N = y(:,param.n+1:end);

clf
% Surface plot
clf
subplot(2,1,1)
surface(t, -param.z,P')
shading flat
xlabel('time (days)')
ylabel('depth (m)')
title('Phytoplankton (mmol N m^-^3)')
colorbar

subplot(2,1,2)
surface(t, -param.z,N')
shading flat
xlabel('time (days)')
ylabel('depth (m)')
title('Nutrients (mmol N m^-^3)')
colorbar
