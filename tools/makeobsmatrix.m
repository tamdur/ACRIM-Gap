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
carrington=0; %1 to use carrington rotation period, 0 to use monthly average
saveString= 'mat_files/obs_23_02_25.mat'; %Name of saved mat file

if PMOD
    paths=dir('observations_pmod/*.txt');
else
    paths=dir('observations/*.txt');
end
paths={paths.name};

stJD = min(dateR); endJD = max(dateR);
stDate = datejd(stJD); endDate = datejd(endJD);
if carrington
    %First make a vector for every carrington period 
    dateM=stDate;
    while juliandate(dateM(end)) < juliandate(endDate)
        dateM = [dateM; dateM(end)+caldays(27)+hours(6)];
    end
    dateMS=dateM(1:end-1);
    dateME=dateM(2:end);
    dateM=mean([dateMS,dateME],2);
else
    %First make a vector for every month in the dataset
    
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
end

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
    [mthVal,mthDate] = dailytomonthly(JD,vals,dateMS,dateME);
    colLabels=[colLabels;string(name)];
    for iB = 1:length(dateMS)
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

%Re-arrange output to be consistent with past output
lI=[5;7;1;2;3;4;6];
oM=oM(:,lI);
colLabels=colLabels(lI);
valM=valM(:,lI);
offsets=offsets(lI);


if ~isempty(saveString)
save(saveString,'oM','dateM','colLabels','valM','offsets');
end
