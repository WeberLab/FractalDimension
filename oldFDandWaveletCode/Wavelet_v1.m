disp('Wavelet Analysis for Fractional Gaussian Noise');
disp(' ');
disp('TYPE ANY KEY to continue to select directory with your data (HEAD and BRIK files)');
disp(' ')
pause
cd(uigetdir)
ls
display(' ')
display('TYPE ANY KEY to continue to select the file you would like to analyze (BRIK file)')
display('note: this may take a while, and MATLAB may continue to display "Pause" even though its actually Busy')
pause
Bbold=uigetfile('*.BRIK', 'Select the BRIK file you would like to analyze');
brainbold = sprintf('%s', Bbold);   % % convert to string
[err,BOLD,Info,ErrMessage]=BrikLoad(brainbold);
BOLD=permute(BOLD,[2,1,3,4]);   % % Rearrange dimensions of array
XX=size(BOLD,1);
YY=size(BOLD,2);
%save('BOLD','BOLD','-v7.3');

BOLD=BOLD(:,:,14:14,:);

ZZ=size(BOLD,3);
TT=size(BOLD,4);

figure(1)
subplot(3,2,1); plot(squeeze(BOLD(32,32,14,:))); title('brain signal') %*
%....
eval(['print -dpsc -r1200 ' num2str(FDdir) '_signal_profile;']); %AMW
close %AMW

cw1 = cwt(squeeze(BOLD(32,32,14,:)),1:9,'gaus2','plot'); % perform wavelet transform of signal using second derivative of the gaussian signal
abcw1=abs(cw1);
% [lmax,I] = localmax(abcw1);
% HEST = wfbmesti(squeeze(BOLD(32,32,14,:)));

% for i=1:9
% for j=1:512
% if lmax(i,j)>0
% lmax(i,j)=500;
% end
% end
% end
% image(lmax)
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
%y(ideb,max_down) = max_down;
%if rowInit<2 , return; end
maxmat=zeros(9,512);
maxmat(1,:)=wmat(1,:);

for jj = ideb+step:step:ifin
    max_curr = find(wmat(jj,:));
    %val_max  = zeros(size(max_curr));
    for k = 1:length(max_up)
        for l=1:length(max_curr)
            if abs(max_curr(l)-max_up(k)) < 2
                maxmat(jj,max_curr(l))=wmat(jj,max_curr(l));
                
            end
        end
    end
    max_up = find(maxmat(jj,:));
end

pl=1;
loc=find(maxmat(1,:));
for pp=1:length(loc)
    max_up=find(maxmat(1,:));
    loca(1,1)=1;
    loca(1,2)=maxmat(1,pp);
    for jj = 2:9
        max_curr = find(maxmat(jj,:));
        for k = 1:length(max_up)
            for l=1:length(max_curr)
                if abs(max_curr(l)-max_up(k)) < 2
                    loca(pl+1,1)=jj;
                    loca(pl+1,2)=max_curr(l);
                    pl=pl+1;
                end
            end
        end
        max_up = find(maxmat(jj,:));
    end
    clear pl
    pl=1;
end

    %         [nul,ind] = min(abs(max_curr-max_up(k)));
    %         val_max(ind) = max_up(k);
    %     end
    %     y(jj,max_curr) = val_max;
    % end
    
    % maxmat=zeros(9,512);
    % zq=zeros(1,3);
    % for go=1:sz
    %     for mgo=1:8
    %         for ngo=1:3
    %             tloc=locs(1,go)-(2-ngo);
    %             if zermatrix(9-mgo, tloc) > 0
    %                 zq(1,ngo)=1;
    %             end
    %         end
    %         if sum(zq) == 1
    %             oneloc=find(zq);
    %         elseif sum(zq) ==2
    
    % maxmat=zeros(9,512);
    % zq=zeros(8,3);
    %
    % for go=1:sz
    %     for ngo=1:3
    %         tloc=locs(1,go)-(2-ngo);
    %         if zermatrix(8, tloc) > 0
    %             zq(8,ngo)=1;
    %         end
    %     end
    %     if sum(zq(8,:)) ==1
    %         loc8=locs(1,go)-(2-find(zq(8,:)));
    %         for ngo=1:3
    %             tloc=loc8-(2-ngo);
    %             if zermatrix(7, tloc) > 0
    %                 zq(7,ngo)=1;
    %             end
    %         end
    %         if sum(zq(7,:))==1
    %             loc7=locs(1,go)-(2-find(zq(7,:)));
    %             for ngo=1:3
    %                 tloc=loc7-(2-ngo);
    %                 if zermatrix(6, tloc) > 0
    %                     zq(6,ngo)=1;
    %                 end
    %             end
    %             if sum(zq(6,:))==1
    %                 loc6=locs(1,go)-(2-find(zq(6,:)));
    %                 for ngo=1:3
    %                     tloc=loc6-(2-ngo);
    %                     if zermatrix(5, tloc) > 0
    %                         zq(5,ngo)=1;
    %                     end
    %                 end
    %                 if sum(zq(5,:))==1
    %                     loc5=locs(1,go)-(2-find(zq(5,:)));
    %                     for ngo=1:3
    %                         tloc=loc5-(2-ngo);
    %                         if zermatrix(4, tloc) > 0
    %                             zq(4,ngo)=1;
    %                         end
    %                     end
    %                     if sum(zq(4,:))==1
    %                         loc4=locs(1,go)-(2-find(zq(4,:)));
    %                         for ngo=1:3
    %                             tloc=loc4-(2-ngo);
    %                             if zermatrix(3, tloc) > 0
    %                                 zq(3,ngo)=1;
    %                             end
    %                         end
    %                         if sum(zq(3,:))==1
    %                             loc3=locs(1,go)-(2-find(zq(3,:)));
    %                             for ngo=1:3
    %                                 tloc=loc3-(2-ngo);
    %                                 if zermatrix(2, tloc) > 0
    %                                     zq(2,ngo)=1;
    %                                 end
    %                             end
    %                             if sum(zq(2,:))==1
    %                                 loc2=locs(1,go)-(2-find(zq(2,:)));
    %                                 for ngo=1:3
    %                                     tloc=loc2-(2-ngo);
    %                                     if zermatrix(1, tloc) > 0
    %                                         zq(1,ngo)=1;
    %                                     end
    %                                 end
    %                                 if sum(zq(1,:))==1
    %                                     loc1=locs(1,go)-(2-find(zq(1,:)));
    %                                     maxy = max([zermatrix(9,locs(1,go)) zermatrix(8,loc8) zermatrix(7,loc7) zermatrix(6,loc6) zermatrix(5,loc5) zermatrix(4,loc4) zermatrix(3,loc3) zermatrix(2,loc2) zermatrix(1,loc1)]);
    %                                     maxmat(9,locs(1,go))=maxy;
    %                                     maxmat(8,loc8)=maxy;
    %                                     maxmat(7,loc7)=maxy;
    %                                     maxmat(6,loc6)=maxy;
    %                                     maxmat(5,loc5)=maxy;
    %                                     maxmat(4,loc4)=maxy;
    %                                     maxmat(3,loc3)=maxy;
    %                                     maxmat(2,loc2)=maxy;
    %                                     maxmat(1,loc1)=maxy;
    %                                 end
    %                             end
    %                         end
    %                     end
    %                 end
    %             end
    %         end
    %     end
    % end
