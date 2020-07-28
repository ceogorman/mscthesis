clearvars

% define EEG data directory
datapath   = 'D:\data\EEG\';
% define directory with helpful extra functions / toolboxes
auxfunpath = 'D:\';

% make fieldtrip toolbox accessible to MATLAB
addpath([auxfunpath,'fieldtrip-20191213\'])
ft_defaults

%% do spectral analysis

% aim: get a baseline alpha power and frequency reading from subs

subj=subject_database('quasirhythm');
% exclude pilot subj
subj(1)=[];

% show feedback plot of alpha power range
quickPlot=false;

for isub=1:numel(subj)
    
    if isub==10 %Missing data
        isub=11
        
    elseif isub==18 %Missing data
        isub=19
        
    else
        
    % load cleaned data
    subjfilename=fullfile(datapath,'cleaned',[subj(isub).id,'_RS_cl.mat']);
    load(subjfilename)
    
    % run spectral analysis
    cfg=[];
    cfg.method='mtmfft'; % multi-taper method for spec smoothing
    % % different outputs
    cfg.output='pow';    % powerspctrm (avg across trials)
    %cfg.output='fourier';% complex sspectrum (by trial)
    % % different tapering, i.e. spectral smoothing methods
    cfg.taper='dpss';    % default tapering method   
    cfg.tapsmofrq=2;     % smoothing width
    %cfg.taper='hann';   % using Hann Window
    %cfg.taper='rectwin';% use rectangular window, i.e. no tapering
    
    cfg.foilim=[1,25];   % freq range of interest
    cfg.pad='nextpow2';  % pad data length to a power of 2
    fr=ft_freqanalysis(cfg,data);
        
    if quickPlot
    figure
    cfg=[];
    cfg.layout  = 'biosemi64.lay';
    cfg.comment = 'no';
    cfg.interactive = 'yes';
    cfg.colormap = 'cool';
    cfg.xlim=[8,12];
    ft_topoplotER(cfg,fr)
    end
    
    % save spectral representation
    outfilename=fullfile(datapath,'fr',[subj(isub).id,'_RS_fr.mat']);
    save(outfilename,'fr');
    end
end



%%
%Plot power spectra for all participants in one figure

subj=subject_database('quasirhythm');
% exclude pilot subj
subj(1)=[];


%Loop gets average power from each participant and adds to new variable
for isub = 1:numel(subj)
    
    if isub==10; %Missing data
        isub=11;
        
    elseif isub==18; %Missing data
        isub=19;
        
    else
    
  subjfilename=fullfile(datapath,'fr',[subj(isub).id,'_RS_fr.mat']);
    load(subjfilename)
    
  
  y = fr.powspctrm;  %y axis is all power spectra - separate lines
  ymod = mean(y,1); %y axis is mean of all rows
  
  sumdat(isub,:) = ymod; %Add data from this loop to master 
  
    end
   
end

x = fr.freq;  %x axis is frequency
sumdat(10,:) = []; %Missing data
sumdat(17,:) = []; %Missing data
powmeans = mean(sumdat,1); %Calculate mean of all participants
transpows = sumdat'; %Rotate rows to columns for plotting
dbtranspows = pow2db(transpows) %Convert to dB
dbpowmeans = pow2db(powmeans) %Convert to dB

figure
p1 = plot(x,dbtranspows,'color',[0.3010 0.7450 0.9330],'LineWidth',0.5) %Plot all 
hold on
xlabel('Frequency (Hz)') %Label x axis
ylabel('Power (dB)') %Label y axis
plot(x,dbpowmeans,'color',[0.8500 0.3250 0.0980],'LineWidth',2) %Add mean to plot
xlim([1 25])
set(gca,'fontsize',22)
hold off

%%
%Individual alpha peaks

subj=subject_database('quasirhythm');
% exclude pilot subj
subj(1)=[];
ploti = 1; %Index for assigning plot locations
plotloc = (1:length(subj)); %Specifies available plot locations


for isub = 1:length(subj)
    
    if isub==10; %Missing data
        %isub=11;
        
    elseif isub==18; %Missing data
        %isub=19;
        
    else

subjfilename=fullfile(datapath,'fr',[subj(isub).id,'_RS_fr.mat']);
    load(subjfilename)
 
 x = fr.freq;  %x axis is frequency
 y = fr.powspctrm;  %y axis is all power spectra - separate lines
 ymod = mean(y,1); %y axis is mean of all rows

%Savitzky-Golay filter
polyOrder = 3;
filterLength = 9;
smSpec = sgolayfilt(-gradient(gradient(ymod)),polyOrder,filterLength); 
 

minalph = 6; %Frequency range to search (Hz)
maxalph = 13;


ymin = find(x==minalph);
ymax = find(x==maxalph);
yrange = ymod(ymin:ymax); %Resize data to only relevant frequencies
xrange = x(ymin:ymax); %Scales x to new y
filtyrange = smSpec(ymin:ymax); %Only relevant frequencies from S-G filter
[pks pklocs] = findpeaks(filtyrange,xrange); %Find peak value and freqs
[pksun pklocsun] = findpeaks(yrange,xrange); %Find unfiltered peaks and freqs

%Take greater amplitude
[v i] = max(pks); %Filtered peaks
pklocs = pklocs(i);
[vun iun] = max(pksun); %Unfiltered peaks
pklocsun = pklocsun(iun);

alphapeaks(isub,1) = pklocs; %Record found peak in Hz


behavsub = subj(isub).id; %Get subject ID as number
behavsub = behavsub(end-1:end);
behavsub = str2num(behavsub);
%Create grid of plots with filtered and unfiltered(/10) spectra for each subject
plottitle = sprintf('Participant %d',behavsub);%Variable for dynamic plot labels
%Variables for filtered and unfiltered peaks
filtpeak = num2str(pklocs);
unfiltpeak = num2str(pklocsun);
if isempty(pklocsun);
    unfiltpeak = 'None';
end

ufylim = max(ymod);%Upper Y limit
ufylim = ufylim/10;
ufylim = ufylim+0.05;
filylim = min(smSpec);%Lower Y limit
filylim = filylim-0.05;
yrange = ufylim-filylim;%Y Range
fpkgloc =ufylim-(yrange*0.1);%Position for filtered peak readout
ufpkgloc = ufylim - (yrange*0.25);%Position for unfiltered peak readout

filplot = subplot(7,6,plotloc(ploti));%Plot filtered data
plot(x,smSpec)
hold on 
%xlim([0 25])
axis([0 25 filylim ufylim])
title(plottitle,'FontSize',8)
%grid


combplot = subplot(7,6,plotloc(ploti));%Overlay unfiltered data
plot(x,(ymod/10));
axis([0 25 filylim ufylim])
legend('Filtered','Unfiltered','Location',[0.68 0.08 0.08 0.08]);
text(20,fpkgloc,filtpeak,'Color',[0 0.4470 0.7410]);
text(20,ufpkgloc,unfiltpeak,'Color',[0.8500 0.3250 0.0980]);
xlabel('Frequency (Hz)','fontsize',8);
ylabel('Power ?','fontsize',8);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',6);
set(gca,'XTickLabelMode','auto');
hold off
set(gca,'fontsize',12)

ploti = ploti+1; %Increment to put next plot in next space
    end
   
end


alphapeaks(10,:) = []; %Missing data
alphapeaks(17,:) = []; %Missing data

%%
%Descriptive statistics on IAFs
meaniaf = mean(alphapeaks)
sdiaf = std(alphapeaks)

histogram(alphapeaks,12)%Distribution of peaks, 12 bins
set(gca,'fontsize',22)
xlabel('IAF (Hz)')
ylabel('N')

