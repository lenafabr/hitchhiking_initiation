function [hit,hitvstime,hitind,contprob,savelinkpos_pf,savepatchpos_pf,wc,lc,inreach,opt] = getcontprob(pos,rotmat,options)
% returns time spent by a Brownian diffuser in contact with a walker
% Inputs: position array, rotation array, options array
% Outputs:
% hit -> number of hits registered
% contprob -> contact probability as a function of time

opt = struct();
% linker length
opt.linklen = 0.1;
% number of linkers
opt.nlinks = 1;
% radius of walker
opt.wrad = 0.1;
% radius of diffuser
opt.drad = 0.1;
% output current step
opt.verbose = 0;
% input linker positions explicitly
opt.linkpos = 0;
% wlc contact probability matrix
opt.wlcprobmat = 0;
% walking indices
opt.walkinds = 1;
% diffusing/tethered target indices
opt.tgtinds = 2;
% sticky endosome
opt.sticky = 0;
% stop when 1st hit registers
opt.stoponhit = 0;
% diffusive linkers
opt.difflinks = 0;
% output linker positions over time?
opt.outputlinkpos = 0;
% use specific receptors on a peroxisome
opt.usepatch = 0;
% input specific receptor positions on the peroxisome
opt.patchpos = 0;
% number of patches?
opt.npatch = 0;
% pex receptor patch radius
opt.patchsize = 0.01;
% diffusive patches?
opt.diffpatch = 0;
% stiff linkers?
opt.stifflink =0;

if(nargin>2); opt = copyStruct(options,opt); end

linklen = opt.linklen;
patchsize = opt.patchsize;
drad = opt.drad;
wrad = opt.wrad;

getwlcprob = 0;
if(any(opt.wlcprobmat(:))&&~opt.stifflink)
	wlcprobmat = opt.wlcprobmat;
	dvalstab = linspace(1e-8,drad+linklen,size(wlcprobmat,1));
	betvalstab = linspace(0,pi,size(wlcprobmat,2));
	getwlcprob = 1;
end

walkinds = opt.walkinds;
nwalk = length(walkinds);

tgtinds = opt.tgtinds;
ntgt = length(tgtinds);

nsnaps = size(pos,1);

if(opt.linkpos)
	%input linker positions
	linkpos_pf = opt.linkpos; % linker positions in particle frame
	nlinks = size(linkpos_pf,1); % number of linkers
else
	%randomly generate linkers
	nlinks = opt.nlinks;
	savelinkpos_pf = zeros(nlinks,4,nwalk,nsnaps);
	linkpos_pf = zeros(nlinks,4,nwalk);
	
	if(opt.difflinks)
		for tc = 1:nsnaps
			for wc = 1:nwalk
				th = acos(2*rand(nlinks,1)-1);
				ph = rand(nlinks,1)*2*pi;
				linkpos_pf(:,:,wc) = [wrad*[sin(th).*cos(ph),sin(th).*sin(ph),cos(th)],1/nlinks*ones(nlinks,1)];
			end
			savelinkpos_pf(:,:,:,tc) = linkpos_pf;
		end
	else
		for wc = 1:nwalk
			th = acos(2*rand(nlinks,1)-1);
			ph = rand(nlinks,1)*2*pi;
			linkpos_pf(:,:,wc) = [wrad*[sin(th).*cos(ph),sin(th).*sin(ph),cos(th)],1/nlinks*ones(nlinks,1)];
		end
		savelinkpos_pf = repmat(linkpos_pf,[1,1,1,nsnaps]);
	end
	
end

usepatch = opt.usepatch;
savepatchpos_pf = 0;
if(usepatch)
	if(any(opt.patchpos(:)))
		% input patch positions
		patchpos_pf = opt.patchpos;
		npatch = size(patchpos_pf,1);
	else
		% randomly generate patches
		npatch = opt.npatch;
		savepatchpos_pf = zeros(npatch,4,ntgt,nsnaps);
		patchpos_pf = zeros(npatch,4,ntgt);
		
		if(opt.diffpatch)
			for tc = 1:nsnaps
				for pc = 1:ntgt
					th = acos(2*rand(npatch,1)-1);
					ph = rand(npatch,1)*2*pi;
					patchpos_pf(:,:,pc) = [drad*[sin(th).*cos(ph),sin(th).*sin(ph),cos(th)],1/npatch*ones(npatch,1)];
				end
				savepatchpos_pf(:,:,:,tc) = patchpos_pf;
			end
		else
			for pc = 1:ntgt
					th = acos(2*rand(npatch,1)-1);
					ph = rand(npatch,1)*2*pi;
					patchpos_pf(:,:,pc) = [drad*[sin(th).*cos(ph),sin(th).*sin(ph),cos(th)],1/npatch*ones(npatch,1)];
			end
			savepatchpos_pf = repmat(patchpos_pf,[1,1,1,nsnaps]);
		end
	end
end

hit = 0;
hitvstime = zeros(nsnaps,1);
inreach = zeros(nsnaps,1);
contprob = zeros(ntgt,1);
hitregistered = 0;
hitind = 0;

for tc = 1:nsnaps
	
	if(~mod(tc,100)&& opt.verbose);fprintf('Processing snapshot %d\n',tc); end
	
	% run over all walkers
	walkcount = 0;
	for wc = walkinds
		walkcount = walkcount+1;
		poswalk = pos(tc,:,wc)'; %column vector
		rotwalk = rotmat(tc,:,:,wc);
				
		% run over all targets (target = peroxisome)
		tgtcount = 0;
		for pc = tgtinds
			tgtcount = tgtcount+1;
			postgt = pos(tc,:,pc)'; % column vector
			rottgt = rotmat(tc,:,:,pc);
			% check if distance betn walker and target is less than minimum distance required for contact
				
			if(norm(poswalk-postgt)<=(linklen+wrad+drad))
				if(opt.sticky) % entire endosome surface is sticky
					hit = 1;
					contprob = 1;
				else % check specific linker tip positions
			
					% run over all endosome linkers
					for lc = 1:nlinks
						if(hitregistered); break; end
						linkpos_lf = poswalk+squeeze(rotwalk)*(savelinkpos_pf(lc,1:3,walkcount,tc)');
						
						if(usepatch)
							patchflag = 0;
							for ptc = 1:npatch
								patchpos_lf = postgt+squeeze(rottgt)*(savepatchpos_pf(ptc,1:3,tgtcount,tc)');
								linkpatchdist = norm(linkpos_lf-patchpos_lf);
								patchflag = (linkpatchdist<(linklen+patchsize));
								if(patchflag==1); break; end
							end
						else
							patchflag = 1;
						end
						
						d = norm(linkpos_lf-postgt);
						minlen = linklen+drad;

						% check if distance betn linker and target is less than minimum reqd
						if(d<=minlen && patchflag)
							
							% get contact probability from the green's function for a wlc (linker)
							if(getwlcprob)
								dotprod = sum((linkpos_lf-poswalk).*(postgt-linkpos_lf))/d/wrad;
								if(abs(dotprod)<1)
									beta = acos(dotprod);
								else
									beta = acos(sign(dotprod));
								end
								dind = round(interp1(dvalstab,1:length(dvalstab),d));
								betind = round(interp1(betvalstab,1:length(betvalstab),beta));
								% wlcdprobmat = tabulated probability of
								% contact
								contprob = wlcprobmat(dind,betind);
							elseif(opt.stifflink)
								linktippos = linkpos_lf+linklen.*(linkpos_lf-poswalk)./norm(linkpos_lf-poswalk);
								if(norm(linktippos-postgt)<=drad)
									contprob = 1;
% 									error('test')
								else
									contprob = 0;
								end
							else
								contprob = 1;
							end
							
							tmp = rand;
							if(tmp<=contprob)
								hit = hit + 1;
								hitvstime(tc) = hit;
% 																error('test')
								hitind = tc;
								
								if(opt.stoponhit)
% 									error('test')
									return
								else
									hitregistered = 1;
								end;
							else
								inreach(tc) = 1;
							end
							
						end
						
					end
				end
				
			end
			
		end
		
	end
	
end

end