%Revised version of makemodernobs.m
%Ted Amdur
%6/15/22

%Develop a set of monthly observations from satellite and proxy
%observations
dateS=getdates;
dateR=dateS.acrimplusfive; %Load range of dates for which to incorporate observations
removeOffsets=1; %1 to remove offsets from observations, 0 to keep native units
satOnly=0; %1 to only use satellites, 0 to use proxies+satellites
PMOD=0;%1 to use Frohlich 2006 PMOD corrections, 0 to use original
saveString= 'mat_files/obs_23_02_01.mat'; %Name of saved mat file

if PMOD
    paths=dir('observations_pmod/*.txt');
else
    paths=dir('observations/*.txt');
end
paths={paths.name};

%First make a vector for every month in the dataset
stJD = min(dateR); endJD = max(dateR);
stDate = datejd(stJD); endDate = datejd(endJD);
stDate=dateshift(stDate,'start','month');
endDate=dateshift(endDate,'end','month');
dateM = stDate;
while juliandate(dateM(end)) < juliandate(endDate)
    dateM = [dateM; dateM(end) + calmonths(1)];
end
dateM = dateM(1:end-1);
dateMS = dateshift(dateM,'start','month');
dateME = dateshift(dateM,'end','month');
dateM=mean([dateMS,dateME],2);

jdS = juliandate(dateMS);
jdE = juliandate(dateME);

%Initialize the output variables
nS=length(paths);
oM = false(length(dateM),nS);
valM = NaN(length(dateM),nS);
colLabels=[];
for ii=1:nS
    tic
    [JD,vals,name] = readobstxt(paths{ii});
    name=string(name);
    %Make monthly data
    [mthVal,mthDate] = dailytomonthly(JD,vals);
    colLabels=[colLabels;string(name)];
    for iB = 1:length(jdS)
        iMth = mthDate >= dateMS(iB) & mthDate < dateME(iB);
        if any(iMth)
            oM(iB,ii)=true;
            valM(iB,ii)=mthVal(iMth);
        end
    end
    toc %display time needed to load each observer
end

if removeOffsets
    for ii = 1:size(valM,2)
        offsets(ii) = nanmean(valM(:,ii));
        valM(:,ii) = valM(:,ii)-offsets(ii);
    end
else
    offsets = [];
end
if satOnly
    for ii=1:size(valM,2)
        if strcmp(colLabels(ii),"BremenMgII")||strcmp(colLabels(ii),"SILSO")
            oM(:,ii)=false;
        end
    end
end

if ~isempty(saveString)
save(saveString,'oM','dateM','colLabels','valM','offsets');
end
