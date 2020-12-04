%%Author: Davide Momi, PhD [momi.davide89@gmail.com],
%%https://twitter.com/davemomi
%%https://davi1990.github.io/


clear all
load('data.mat');
d=abs(data_exp);

%%
samps1 = 1;
samps2 = 499;

%bootstrap:
for bb=1:500 %parameters.bootstrap % number of bootstrap (e.g. 500)
    bb
    randtrials=rand([1,size(d,3)]);%vector of length (number of trials) containing random real numbers between 0 and 1
    randonsamp=samps1+ round(randtrials.*(length(samps1:1:samps2-1)-1));% vector of length (number of trials) containing random integer numbers between samp1 and samp2
    for tt=1:1:length(randonsamp)
        H1(:,tt)=d(:,randonsamp(tt),tt);%matrix with dimensions (channels x 1 x trials)
    end
    DataFile.bootstrap(:,bb)=mean(H1,2); %matrix with dimensions (channels x number of bootstrap)
end


avgchans=1;

GH= sqrt( sum( DataFile.bootstrap (avgchans,:).^2,1)/length(avgchans) );%surrogate gmfp

sGH=sort(GH);
parameters.alpha = 3.3333e-05 %p-value threshold (e.g. 0.01)
boot1envelope=repmat(sGH(round(size(GH,2)*(1-parameters.alpha))),1,size(d,2));
clear sGH

times = -500:999

figure;
plot(times,d,'-k')
hold on
xa=plot(times,boot1envelope,'--r');





%%
npnts = 500;
time  = linspace(1,500,npnts)'; % force column vector
nPerms  = 1000; % number of iterations in permutation testing
pval    =  .05; % p-value threshold

sigThresh = norminv(1-pval/2); % note: two-tailed!


data1= d(:,1:500);
data2 = d(:,501:1000);


data3d = cat(2,data1,data2);


% generate true condition labels
condlabels = (1:npnts*2)>npnts;


% initialize
permutedDiffs = zeros(npnts,nPerms);

for permi = 1:nPerms
    
    % shuffle condition label vector
    fakeconds = condlabels( randperm(2*npnts) );
    
    % compute and store difference time series
    mean1 = mean( data3d(:,fakeconds==0),1 );
    mean2 = mean( data3d(:,fakeconds==1),1 );
    permutedDiffs(:,permi) = mean2-mean1;
end




%%
% initialize cluster sizes from permutation
nPerms = 1000;
clustsizes = zeros(nPerms,1);

for permi=1:nPerms
    
    % compute z-score difference
    zdiffFake = (permutedDiffs(:,permi)-mean(permutedDiffs,2)) ./ std(permutedDiffs,[],2);
    
    % threshold
    zdiffFake( abs(zdiffFake)<sigThresh ) = 0;
    
    % identify clusters
    islands = bwconncomp( logical(zdiffFake) );
    
    % find cluster sizes
    clustNs = cellfun(@length,islands.PixelIdxList);
        if isempty(clustNs);
            clustNs=0;
        end
            
    clustsizes(permi) = max(clustNs);
end

% compute cluster threshold
clustthresh = prctile(clustsizes,100-pval*100);


% show distribution of cluster sizes
figure(2), clf
histogram(clustsizes)
hold on
plot([1 1]*clustthresh,get(gca,'ylim'),'r--','linew',3)
xlabel('Cluster size (time points)'), ylabel('Count')


%
boot1envelope=boot1envelope(1:500);
times = times(501:1000);

%


mean1 = mean( data3d(:,condlabels==0),1 );
mean2 = mean( data3d(:,condlabels==1),1 );
real_Diff = mean2-mean1;
% 
media = squeeze(mean(mean(permutedDiffs,2)));
standard = squeeze(mean(std(permutedDiffs',1)));
% 
zdiff  =abs((real_Diff - media') ./ (standard));



% recompute thresholded time series
zthresh = zdiff;
zthresh( abs(zthresh)<boot1envelope(1) ) = 0;



figure(3), clf
subplot(311), hold on
plot(times,zdiff,'-k')
hold on
xa=plot(times,boot1envelope,'--r');
title('Statistical results, uncorrected')

% find islands
islands = bwconncomp( logical(zthresh) );

% find and remove any subthreshold islands
for ii=1:islands.NumObjects
    if numel(islands.PixelIdxList{ii})<clustthresh
        zthresh(islands.PixelIdxList{ii}) = 0;
    end
end

% now plot that
subplot(312), hold on
plot(times,zdiff,'k')
plot(times(logical(zthresh)), zthresh(logical(zthresh)),'yo','markerfacecolor','r')
xa=plot(times,boot1envelope,'--r');
ylabel('Z value')
title('Statistical results, cluster-corrected')
legend({'Difference';[ 'p<' num2str(pval) ', 2-tailed' ]})




% now plot that nicer
subplot(313), cla, hold on
plot(times,zdiff,'k')
xa=plot(times,boot1envelope,'--r');

islands = bwconncomp( logical(zthresh) );
for ii=1:islands.NumObjects
    iidx = islands.PixelIdxList{ii};
    patch([ times(iidx)' times(iidx(end:-1:1))' ],[zdiff(iidx)' zeros(1,numel(iidx))'],'r','edgecolor','none');
end

xlabel('Time (s)'), ylabel('Abs Source Amplitute')
title('Statistical results, corrected')
legend({'Difference';[ 'p<' num2str(pval) ', 2-tailed' ]})
