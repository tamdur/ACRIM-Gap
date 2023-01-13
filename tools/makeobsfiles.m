%Generate observational files for each datasource, both satellite and
%observational
% 
% Ted Amdur
% 6/3/22
clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load satellite files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dateS=getdates;
dateR=dateS.all;
[compositeObs,ds] = getdailyobs(datejd(min(min(dateR))),datejd(max(max(dateR))));
sats = unique(ds.ID);
%colLabels = [sats2';"SILSO";"BremenMgII"];
makeinstruments
instruments=sat; %Hold the instrument structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add proxy instrument files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
satL=length(instruments);
%sunspots
instruments(satL+1).ID ="SILSO";
instruments(satL+1).version="v2";
vInd=~isnan(compositeObs.spotsInt);
instruments(satL+1).JD=compositeObs.day(vInd)+0.5;
instruments(satL+1).TSI=compositeObs.spotsInt(vInd);

%Mg-II
instruments(satL+2).ID ="BremenMgII";
instruments(satL+2).version="v5";
vInd=~isnan(compositeObs.mg);
instruments(satL+2).JD=compositeObs.day(vInd)+0.5;
instruments(satL+2).TSI=compositeObs.mg(vInd);

for ii=1:length(instruments)
    colLabels(ii)=string(strrep(char(instruments(ii).ID),'/','_'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Write text files for each 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:length(instruments)
    filePath=[char(colLabels(ii)) '.txt'];
    fid=fopen(filePath,'w');
    
    %Print a line for the title
    fprintf(fid,'%s\n',char(instruments(ii).ID));
    
    %Print a line for the version
    vLine=['Version: ' char(instruments(ii).version)];
    fprintf(fid,'%s\n',vLine);
    
    %Print a line for when created
    dLine=['Created: ' char(datetime('today'))];
    fprintf(fid,'%s\n',dLine);
    
    %Print output
    JD=instruments(ii).JD;Value=instruments(ii).TSI;
    X=[JD Value];
    fprintf(fid,'%s %12s \n','JD','Value');
    fprintf(fid,'%7.1f %f\n',X');
    fclose(fid);
end



function timeMatchedVals = timematch(t,tvals,vals)
[~,uI] = unique(tvals); vals = vals(uI); %Make sure there's no repeats
tvals = floor(tvals(uI)); %Floor in case JD ends on 0.5
valI = (1:length(vals))'; %index of values
leftInd = (tvals-t(1)+1)'; %Index of insertion
valid = leftInd > 0 & leftInd <= max(t); %vals that can be added to timeMatchedVals
rightInd = valI(valid); %Select the values to be placed into timeMatchedVals
leftInd = leftInd(valid); %Select where to place timeMatchedVals
timeMatchedVals = NaN(size(t)); %Create timeMatchedVals vector of size t
timeMatchedVals(leftInd) = vals(rightInd); %Insert
end


function obsMG = getmg
%import file as table
importmgiiupdated
obsMG = MgIIcomposite211129;
obsMG.datetime = datetime(floor(obsMG.fracYr),obsMG.Month,obsMG.Day);
obsMG.JD(:) = juliandate(obsMG.datetime);
obsMG.MG(:) = obsMG.mgii;
end
