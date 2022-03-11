function dydt = PmodelDerivNPSeason(t,y,param)

P = y(1:param.n);
N = y((param.n+1):end);

%--------------------PHYTOPLANKTON-----------------------------------------
% Interior fluxes:
for i = 2:param.n
    Ja(i) = param.u * P(i-1); % advective flux is the velocity x conc. in previous cell
    Jd(i) = -param.D * (P(i)-P(i-1))/param.dz; % diffusive fluxes
end

% Boundary fluxes:
Ja(1) = 0; % no input from surface
Jd(1) = 0; % surface
Ja(param.n+1) = 0; % no flux into the bottom
Jd(param.n+1) = 0; % no flux at bottom

% Total flux:
J = Ja + Jd;

%----------------------NUTRIENTS-------------------------------------------
% Interior fluxes:
for i = 2:param.n
    Jan(i)= 0; %advective fluxes (same as P)
    Jdn(i) = -param.D*(N(i)-N(i-1))/param.dz; %diffusive fluxes (same as P)
end

% Boundary fluxes:
Jan(1) = 0; % no input from surface
Jdn(1) = 0; % surface
Jan(param.n+1) = 0; % bottom
Jdn(param.n+1) = - param.D * ((param.nz - N(end))/param.dz); % difussive flux at bottom

% Total flux:
Jn = Jan + Jdn;

%--------------------------N AND P DYNAMICS--------------------------------
% Phytoplankton growth:
I = CalcLightSeason(P', param,t);

% Liebig's law phytoplankton growth rate:
mu = param.gmax*min(N./(param.H_N + N),(I./(param.H_I + I)));

% Updating growth equation:
for i = 1:param.n
    dPdt(i) = -((J(i+1)-J(i))/param.dz) + mu(i).*P(i) - param.m.*P(i);
end

% Updating nutrient equation:
for i = 1:param.n
    dNdt(i) = -((Jn(i+1)-Jn(i))/param.dz) - param.alpha*mu(i).*P(i) + param.eps*param.alpha*param.m.*P(i);
end

% Make dPdt and dNdt column vectors:
dydt = [dPdt dNdt]';

end
