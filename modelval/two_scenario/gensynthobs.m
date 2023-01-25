function valS = gensynthobs(scenario,dateM,Ainit,epsilon,rho)
%GENSYNTHOBS Generate a set of synthetic observations for ACRIM-Gap BTSI
%testing.
% INPUTS:
%           scenario: char vector, either 'ACRIM' or 'PMOD', that informs
%           which baseline TSI scenario to use. Error otherwise
%           Ainit: assumed true observation matrix, from
%           initobservationmodelparams.m
%           epsilon: assumed instrument noise, from initobservationmodelparams.m
%           rho: assumed autocorrelation parameter for satellite noise,
%           from initobservationmodelparams.m
% OUTPUTS:
%           valS: valM observation array with synthetic data generated
%           using input assumptions. Proxies include offset and synthetic Gaussian
%           noise, and satellites include offset, linear drift, and AR(1)
%           autocorrelated residual error

%Load TSI series used for the assumed truth
tsi = twotsiseries;

if strcmp(scenario,'ACRIM')
    tsi=tsi.ACIRM;
elseif strcmp(scenario,'PMOD')
    tsi=tsi.PMOD;
else
    error('Enter "ACRIM" or "PMOD" for the scenario input')
end
datetsi=tsi.date;
    

end

