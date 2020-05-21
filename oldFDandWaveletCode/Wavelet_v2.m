% disp('Wavelet Analysis for Fractional Gaussian Noise');
% disp(' ');
% disp('TYPE ANY KEY to continue to select directory with your data (HEAD and BRIK files)');
% disp(' ')
% pause
% cd(uigetdir)
% ls
% display(' ')
% display('TYPE ANY KEY to continue to select the file you would like to analyze (BRIK file)')
% display('note: this may take a while, and MATLAB may continue to display "Pause" even though its actually Busy')
% pause
% Bbold=uigetfile('*.BRIK', 'Select the BRIK file you would like to analyze');
% brainbold = sprintf('%s', Bbold);   % % convert to string
% [err,BOLD,Info,ErrMessage]=BrikLoad(brainbold);
% BOLD=permute(BOLD,[2,1,3,4]);   % % Rearrange dimensions of array
% XX=size(BOLD,1);
% YY=size(BOLD,2);
% save('BOLD','BOLD','-v7.3');
% 
% ZZ=size(BOLD,3);
% TT=size(BOLD,4);

%cw1 = cwt(squeeze(BOLD(32,32,14,:)),1:9,'gaus2','plot'); % perform wavelet transform of signal using second derivative of the gaussian signal
%cw1 = cwt((FBM),1:9,'gaus2','plot');

cw1 = cwt(sampsig,1:9,'gaus2');
abcw1=abs(cw1);

wmat=zeros(9,512);
for mgo=1:9
    [peaks,locs]=findpeaks(abcw1(mgo,:));
    sz=size(peaks);
    sz=sz(1,2);
    for go=1:sz
        wmat(mgo,locs(1,go))=peaks(1,go);
    end
end

ideb = 1 ; step = 1; ifin = 9;
max_up=find(wmat(ideb,:));

maxmat=zeros(9,512);
maxmat(1,:)=wmat(1,:);

for jj = ideb+step:step:ifin
    max_curr = find(wmat(jj,:));
    
    for k = 1:length(max_up)
        for l=1:length(max_curr)
            if abs(max_curr(l)-max_up(k)) < 2
                maxmat(jj,max_curr(l))=wmat(jj,max_curr(l));
            end
        end
    end
    max_up = find(maxmat(jj,:));
end

wmat=zeros(9,512);
loc=find(maxmat(1,:));
for pp=1:length(loc)
    indeces=[1 loc(pp) maxmat(1,loc(pp))];
    for jj=2:9
        if size(indeces,1) >= jj-1
            up=indeces(jj-1,2);
            curr=find(maxmat(jj,:)); %
            for l=1:length(curr)
                if abs(curr(l)-up) == 0
                    newind=[jj curr(l) maxmat(jj,curr(l))];
                    indeces=[indeces; newind];
                elseif abs(curr(l)-up) == 1
                    newind=[jj curr(l) maxmat(jj,curr(l))];
                    indeces=[indeces; newind];
                end
            end
        end
    end
    newmax=max(indeces(:,3));
    rows=size(indeces,1);
    for ii=1:rows
        wmat(indeces(ii,1),indeces(ii,2))=newmax;
    end
    clear indeces
    clear newmax
end

z=zeros(9,2);
for pp=1:9
    z(pp,1)=log(sum(wmat(pp,:)));
    z(pp,2)=log(pp);
end

brob=robustfit(z(:,2),z(:,1));
hurst=1-brob(2);


figure(1)
subplot(6,1,1); plot(sampsig); title('brain signal')
subplot(6,1,2); imagesc(flipud(cw1)); title('gaus2 wavelet')
subplot(6,1,3); imagesc(flipud(abcw1)); title('|gaus2 wavelet|')
subplot(6,1,4); imagesc(flipud(maxmat)); title('local max')
subplot(6,1,5); imagesc(flipud(wmat)); title('supremum')
subplot(6,1,6); scatter(z(:,2),z(:,1),'filled'); grid on; hold on
plot(z(:,2),brob(1)+brob(2)*z(:,2),'g','LineWidth',2); title('log sum supremum v log scale')
%eval(['print -dpsc -r1200 ' num2str(FDdir) '_signal_profile;']);
%close %AMW