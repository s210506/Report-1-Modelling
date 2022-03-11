function I = CalcLightSeason(P, param,t)
% Light

% Integral to account for P dampning
int = cumsum(P.*param.dz)*param.k;

% Seasonal forcing:
season = 1-0.8 * cos(2*pi*(t/365));

% Light function:
%I = param.I0*exp(-param.Kbg*param.z-int)*(1-0.8 * cos(2*pi*(t/365)));
I = param.I0*exp(-param.Kbg*param.z-int)*season;

% Making light a column vector:
I = I';
end

