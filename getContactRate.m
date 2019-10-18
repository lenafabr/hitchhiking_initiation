%% use simulation snapshots to calculate rate of initial contact with a linker tip
% for a single value of parameters, but many trajectory iterations
% this script calculates linker contact rates for varying linker lengths or flexibility
%% load in tabulated contact probabilities
clear;clc;

% ---------------
% Code for generating distribution function for wormlike chains and linker  contact probability matrices can be provided upon request (to smogre@ucsd.edu)
% ----------------

% uncomment following line to load tabulated contact probabilities for varying length
% the array linklens contains the linker lengths for which the contact probability is calculated
load('./wlcfiles/contactprobmat_Rpex0pt1_lp0pt1_ldvals20nmto200nm.mat')

% uncomment following line to load tabulated contact probabilities for varying flexibility
% load('./wlcfiles/contactprobmat_Rpex0pt1_lp10nmto1um_ld100nm.mat'); 

%% load in simulation trajectory data
rhoc = 3; %number of carrier organelles
trkc = 5; %number of microtubules
tethstr = 'diff'; %tethering state
load(sprintf('./posfile_%dendo_%dMT_%s.mat',rhoc,trkc,tethstr));
%% set up options
nlinks = 5; %number of linkers
linklen = 0.1; %linker length
lp = 0.1; %linker flexibility

options = struct();
options.stoponhit = 1;
options.difflinks = 0; %do the linkers diffuse on the carrier surface?
options.drad = drad; %peroxisome radius for calculating linker contact
options.nlinks = nlinks; 
options.walkinds = walkinds; %indices of walking particles
options.tgtinds = npart; %number of particles

% uncomment following lines to vary linker length
[~,lenc] = min(abs(linklens-linklen));
options.linklen = linklens(lenc);
options.wlcprobmat = contprobmat{lenc};

% uncomment following lines to vary linker flexibility
% [~,lpc] = min(abs(lpvals-lp));
% options.linklen = linklen;
% options.wlcprobmat = contprobmat{lpc};

%% get initial contact times for all the iterations
newtmp = zeros(size(pos,1),1);
parfor tc = 1:size(pos,4) % cycle over trials
    % tmp contains 1 at the first timepoint where a hit occured
    % zeros elsewhere
    [~,tmp] = linkerContactSim(pos(:,:,:,tc),rotmat(:,:,:,:,tc),options);
    newtmp = newtmp+tmp;
end
hitvtime = newtmp/size(pos,4);
	
%% plot cumulative fraction hit over time
cumhit = cumsum(hitvtime);

%fit to double exponential
func = @(c,t) 1 - c(1)*exp(-t/c(2))-c(3)*exp(-t/c(4));

lb = [0,0.01,0,0.01];
ub = [1,100,1,100];
cfit = lsqcurvefit(func,[0.5,0.1,0.5,1],tvals,cumhit,lb,ub);

plot(tvals,cumhit,tvals,func(cfit,tvals))

khit = 1/(cfit(1)*cfit(2)+cfit(3)*cfit(4))