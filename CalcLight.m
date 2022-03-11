function I = CalcLight(P, param)
% Light

% Integral to account for P dampning
int = cumsum(P.*param.dz)*param.k;

% Light function:
I = param.I0*exp(-param.Kbg*param.z-int);

% Making light a column vector:
I = I';
end

