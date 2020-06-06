function sralt_plot(destDir, numImage, canonicalImageSize, layout)

% initial input images
load(fullfile(destDir, 'original-sralt.mat'), 'D') ;

% alignment results
load(fullfile(destDir, 'final-sralt.mat'), 'Do','A','E') ;

%% display
% layout
if nargin < 4
    xI = ceil(sqrt(numImage)) ;
    yI = ceil(numImage/xI) ;

    gap = 2;
    gap2 = 1; 
else
    xI = layout.xI ;
    yI = layout.yI ;

    gap = layout.gap ;
    gap2 = layout.gap2 ; 
end
container = ones(canonicalImageSize(1)+gap, canonicalImageSize(2)+gap); 
% white edges
bigpic = cell(xI,yI); 

% D 
for i = 1:xI
    for j = 1:yI
        if yI*(i-1)+j > numImage
            bigpic{i,j} = ones(canonicalImageSize(1)+gap, canonicalImageSize(2)+gap);
        else
            container ((gap2+1):(end-gap2), (gap2+1):(end-gap2)) = reshape(D(:,yI*(i-1)+j), canonicalImageSize);
            bigpic{i,j} = container;
        end
    end
end
figure(1)
imshow(cell2mat(bigpic),[],'DisplayRange',[0 max(max(D))],'Border','tight')
title('Input images') ;
savefile = [destDir  '\original-sralt.fig'];
saveas(1,savefile,'fig');
            

% Do¡¡
for i = 1:xI
    for j = 1:yI
        if yI*(i-1)+j > numImage
            bigpic{i,j} = ones(canonicalImageSize(1)+gap, canonicalImageSize(2)+gap);
        else
            container ((gap2+1):(end-gap2), (gap2+1):(end-gap2)) = reshape(Do(:,yI*(i-1)+j), canonicalImageSize);
            bigpic{i,j} = container;
        end
    end
end
figure(2)
imshow(cell2mat(bigpic),[],'DisplayRange',[0 max(max(Do))],'Border','tight')
title('Aligned images') ;
savefile = [destDir  '\Do-sralt.fig'];
saveas(2,savefile,'fig');


% A¡¡
for i = 1:xI
    for j = 1:yI
        if yI*(i-1)+j > numImage
            bigpic{i,j} = ones(canonicalImageSize(1)+gap, canonicalImageSize(2)+gap);
        else
            container ((gap2+1):(end-gap2), (gap2+1):(end-gap2)) = reshape(A(:,yI*(i-1)+j), canonicalImageSize);
            bigpic{i,j} = container;
        end
    end
end
figure(3)
imshow(cell2mat(bigpic),[],'DisplayRange',[0 max(max(A))],'Border','tight')
title('Aligned images adjusted for sparse errors with L') ;
savefile = [destDir  '\L-sralt.fig'];
saveas(3,savefile,'fig');

% E¡¡
for i = 1:xI
    for j = 1:yI
        if yI*(i-1)+j > numImage
            bigpic{i,j} = ones(canonicalImageSize(1)+gap, canonicalImageSize(2)+gap);
        else
            container ((gap2+1):(end-gap2), (gap2+1):(end-gap2)) = reshape(E(:,yI*(i-1)+j), canonicalImageSize);
            bigpic{i,j} = container;
        end
    end
end
figure(4)
imshow(abs(cell2mat(bigpic)),[],'DisplayRange',[0 max(max(abs(E)))],'Border','tight')
title('Sparse corruptions in the aligned images') ;
savefile = [destDir  '\E-sralt.fig'];
saveas(4,savefile,'fig');

figure(5)
subplot(1,3,1)
imshow(reshape(sum(D,2), canonicalImageSize),[])
title('average of unaligned D')
subplot(1,3,2)
imshow(reshape(sum(Do,2), canonicalImageSize),[])
title('average of aligned D')
subplot(1,3,3)
imshow(reshape(sum(A,2), canonicalImageSize),[])
title('average of A')
savefile = [destDir  '\average-sralt.fig'];
saveas(5,savefile,'fig');