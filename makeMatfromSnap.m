%% load sim data and save to file
clear;clc;
rhoc = 3; %number of carrier organelles
trkc = 5; %number of microtubules
tethstr = 'diff'; %tethering state
savefilestr = sprintf('./posfile_%dendo_%dMT_%s.mat',rhoc,trkc,tethstr);
snapfilestr = sprintf('./BDsim_%dendo_%dMT_%s.snap.out',rhoc,trkc,tethstr);

if(exist(snapfilestr,'file'))
	fprintf('Extracting rhoc = %d, trkc = %d, %s\n',rhoc,trkc,tethstr)
	[pos,rotmat,nsnap,npart,ntrials,tvals,domrad,domlen,walkinds,tethinds,sprpos,pinfo]...
		= loadSnapshot(snapfilestr);
	sterrad = pinfo(3);
	nteth = length(tethinds);
	nwalk = length(walkinds);
	save(savefilestr,'-v7.3')
else
	fprintf('rhoc = %d, trkc = %d, %s does not have a snap file\n',rhoc,trkc,tethstr)
end