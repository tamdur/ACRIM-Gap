function [compositeObs,ds] = getdailyobs(timeStart,timeEnd)

%Assimilate all the observational data needed for the gibbs sampling model. Heavily
%drawn from compsoitespot_coordinate_update.m and related scripts.

timeStart = floor(juliandate(timeStart));
timeEnd = floor(juliandate(timeEnd));

%Create a JD vector that indexes all the days that we'll use
compositeObs.day = (timeStart:timeEnd)';
compositeObs.ind = (1:length(compositeObs.day))';

%Load PMOD TSI record
importPMOD
%note that 1/1/1980 corresponds to 1.5 in epoch1980
PMOD.JD = PMOD.epoch1980 + juliandate(datetime(1980,01,01))-1.0; 
%Make NaNs for invalid observations
PMODNew = PMOD.TSI_new_;
PMODOld = PMOD.TSI_old_;
PMODNew(PMOD.TSI_new_ < 0) = NaN;
PMODOld(PMOD.TSI_old_ < 0) = NaN;
PMOD.JD = floor(PMOD.JD);
[~,iiPMOD] = unique(PMOD.JD);
PMOD = PMOD(iiPMOD,:);
PMODNew = PMODNew(iiPMOD);

%Load mgII
obsMG = getmg;

%load radio
load('obsRadioUpdated.mat') %Comes from penticton record
obsRadio.JD = floor(obsRadio.JD);
    
%load int sunspot number
makesilso
obsSilso.JD = floor(obsSilso.JD);

%cut down the sunspot data to days with corresponding TSI values
compositeObs.mg = timematch(compositeObs.day,obsMG.JD,obsMG.MG);
compositeObs.mgError = timematch(compositeObs.day,obsMG.JD,obsMG.uncertainty);
compositeObs.radio = timematch(compositeObs.day,obsRadio.JD, obsRadio.radio);
compositeObs.spotsInt = timematch(compositeObs.day,obsSilso.JD,obsSilso.spotnum);
compositeObs.PMODNew = timematch(compositeObs.day,PMOD.JD,PMODNew);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load satellites
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ds = makedailysatdataset(timeStart,timeEnd);
end

function obsMG = getmg
%import file as table
importmgiiupdated
obsMG = MgIIcomposite;
obsMG.datetime = datetime(floor(obsMG.fracYr),obsMG.Month,obsMG.Day);
obsMG.JD(:) = juliandate(obsMG.datetime);
obsMG.MG(:) = obsMG.mgii;
end

function ds = makedailysatdataset(startDate,endDate)
%makedailysatdataset Produce dataset object, where each row is a TSI
%observation with corresponding metadata and solar conditions for that day.
%   makelmedataset produces a dataset array for all daily observations from
%   the 10 satellites used in the instrumental record of makeInstruments.m
%   and the sponspot records compiled in compositespot_coordinate_update. 
%
%           
%   OUTPUTS: 
%           ds: dataset for which each row is a satellite observation, each
%           column is an associated variable
%    
makeinstruments
JD = [];
TSI = [];
ID = [];
satAge = [];
for ii = 1:length(sat)
%     if strcmp(sat(ii).ID, 'NIMBUS/HF') %start Nimbus/HF in 1982 after issues
%         dAvail = find(sat(ii).JD<= endDate ...
%         & sat(ii).JD >= juliandate(datetime(1982,01,01)));
%     else
        dAvail = find(sat(ii).JD<= endDate ...
        & sat(ii).JD >= startDate); %keep within compositeObs
%     end
    noNan = ~isnan(sat(ii).TSI(dAvail)); %remove missing vals
    dAvail = dAvail(noNan);
    %Only accept one value per day
    [~,iiJD] = unique(sat(ii).JD(dAvail));
    dAvail = dAvail(iiJD);
    TSI = [TSI; sat(ii).TSI(dAvail)];
    JD = [JD; sat(ii).JD(dAvail)];
    thisSat = repmat(string(sat(ii).ID), length(dAvail),1);
    ID = [ID; cellstr(thisSat)];
    if  ~isempty(dAvail)
        satAge = [satAge; yearfraction(datejd(sat(ii).JD(dAvail)),datejd(sat(ii).JD(dAvail(1))))];
        %satAge = [satAge; yearfraction(juliandate(sat(ii).JD(dAvail)),juliandate(sat(ii).JD(dAvail(1))))];
    end
end

% Assemble the vectors into a dataset object
yearf = yearfraction(datejd(JD),datejd(startDate));
%yearf = yearfraction(juliandate(JD),juliandate(startDate));
year = floor(yearf);
% TSImin1 = TSI; TSImin1(2:end) = TSImin1(1:end-1);
instruments = dataset({ID,'ID'},...
    {yearf,'YEARF'},...
    {JD,'JD'},...
    {year,'YEAR'},...
    {satAge,'satAge'},...
    {TSI,'TSI'});
ds = instruments;
ds.ID = nominal(ds.ID);
ds.YEAR = nominal(ds.YEAR);
end

function timeMatchedVals = timematch(t,tvals,vals)
[~,uI] = unique(tvals); vals = vals(uI); %Make sure there's no repeats
tvals = floor(tvals(uI)); %Floor in case JD ends on 0.5
valI = (1:length(vals))'; %index of values
leftInd = (tvals-t(1)+1)'; %Index of insertion
valid = leftInd > 0 & leftInd <= length(t); %vals that can be added to timeMatchedVals
rightInd = valI(valid); %Select the values to be placed into timeMatchedVals
leftInd = leftInd(valid); %Select where to place timeMatchedVals
timeMatchedVals = NaN(size(t)); %Create timeMatchedVals vector of size t
timeMatchedVals(leftInd) = vals(rightInd); %Insert
end




