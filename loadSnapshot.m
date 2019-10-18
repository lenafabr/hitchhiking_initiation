function [pos,rotmat,nsnap,npart,ntrials,tvals,domrad,domlen,walkinds,tethinds,sprpos,pinfo] = loadSnapshot(filename)

% load data on particle positions and rotations from a snapshot file
% ------------
% outputs:
% ------------
% pos(i,p,s,t) has the i coordinate of the p particle in the s snapshot in simulation trial t
% rotmat(i,j,p,s) has the i,j matrix component of the orientation matrix
% for the p particle in the s snapshot in simulation trial t
% nsnap is number of snapshots
% npart is number of particles
% ntrials is the bumber of trials
% tvals are the times of the snapshots
% domrad is the domain radius
% domlen is the domain half-length
% walkinds are indices p of the walking particles
% tethinds are indices p of tethered particles (if any)
% sprpos gives the spring position if particles are walking
% pinfo is generic info about particles (fric, fricr, srad, crad)
% (currently assumed constant for all particles and all time)

data = dlmread(filename);

ntrials = data(1,1);
npart = data(1,2);
nwalk = data(1,3);
nteth = data(1,4);
pinfo = data(1,6:end);
if(nwalk == 0); walkinds = []; else; walkinds = pinfo(4+(1:nwalk)); end
if(nteth == 0); tethinds = []; else; tethinds = pinfo(4+nwalk+(1:nteth)); end
domrad = pinfo(end-4);
domlen = pinfo(end-3);

nsnap = size(data,1)/(ntrials*(npart+1));
np1 = npart+1;
pos = zeros(nsnap,3,npart,ntrials);
sprpos = zeros(nsnap,3,nwalk+nteth,ntrials);
rotmat = zeros(nsnap,3,3,npart,ntrials);
tvals = zeros(nsnap,1);

for tc = 1:ntrials
  nrows = np1*nsnap;
  tmpdata = data((tc-1)*nrows+(1:nrows),:);
  for cc = 1:nsnap
    
    snapshot = tmpdata(np1*(cc-1)+2:np1*cc,:);
    tvals(cc) = tmpdata(np1*(cc-1)+1,5);

    pos(cc,:,:,tc) = snapshot(:,2:4)';
    for pc = 1:npart
      rotmat(cc,:,:,pc,tc) = reshape(snapshot(pc,5:13),3,3);
    end
    if(nwalk>0 || nteth>0)
      sprpos(cc,:,:,tc) = snapshot([walkinds,tethinds],14:16)';
    end
  end
end

end