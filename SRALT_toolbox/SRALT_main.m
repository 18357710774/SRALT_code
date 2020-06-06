function [Do, A, E, xi, numIterOuter, numIterInner] = SRALT_main(fileNames, transformations, numImages, para, destDir)

if ~isfield(para,'inner_tol');       para.inner_tol = 1e-7;      end
if ~isfield(para,'inner_maxIter');       para.inner_maxIter = 1000;      end
if ~isfield(para,'muBound')     
    if method == 0
        para.muBound = 1e+30; 
    else
        para.muBound = 1e-30; 
    end
end
if ~isfield(para,'rho0');       para.rho0 = 1.25;      end


%% read and store full images
fixGammaType = 1 ;

if ~fixGammaType
    if exist(fullfile(rootPath, userName, 'gamma_is_ntsc'), 'file')
        gammaType = 'ntsc' ;
    elseif exist(fullfile(rootPath, userName, 'gamma_is_srgb'), 'file')
        gammaType = 'srgb' ;
    elseif exist(fullfile(rootPath, userName, 'gamma_is_linear'), 'file')
        gammaType = 'linear' ;
    else
        error('Gamma type not specified for training database!  Please create a file of the form gamma_is_*') ;
    end
else
    gammaType = 'linear' ;
end

sigma0 = 2/5 ;
sigmaS = 1 ;

deGammaTraining = true ;

I0 = cell(para.numScales,numImages) ; 
I0x = cell(para.numScales,numImages) ; 
I0y = cell(para.numScales,numImages) ;

for fileIndex = 1 : numImages
    
    currentImage = double(imread(fileNames{fileIndex}));
    if size(currentImage,3) > 1,   currentImage = currentImage(:,:,2);            end
    
    if deGammaTraining,      currentImage = gamma_decompress(currentImage, gammaType); end

    currentImagePyramid = gauss_pyramid( currentImage, para.numScales,...
        sqrt(det(transformations{fileIndex}(1:2,1:2)))*sigma0, sigmaS );
        
    for scaleIndex = para.numScales:-1:1
        I0{scaleIndex,fileIndex} = currentImagePyramid{scaleIndex};
        I0_smooth = I0{scaleIndex,fileIndex};
        I0x{scaleIndex,fileIndex} = imfilter( I0_smooth, (-fspecial('sobel')') / 8 );
        I0y{scaleIndex,fileIndex} = imfilter( I0_smooth,  -fspecial('sobel')   / 8 );
    end   
end

%% get the initial input images in canonical frame
imgSize = para.canonicalImageSize ;   
xi_initial = cell(1,numImages) ;
for i = 1 : numImages
    if size(transformations{i},1) < 3
        transformations{i} = [transformations{i} ; 0 0 1] ;
    end
    % xi_initial：[\tau_1,...,\tau_n];
    xi_initial{i} = projective_matrix_to_parameters(para.transformType,transformations{i});  
end

D = zeros(imgSize(1)*imgSize(2), numImages); 
for fileIndex = 1 : numImages
    % transformed image
    Tfm = fliptform(maketform('projective',transformations{fileIndex}'));  
    I   = vec(imtransform(I0{1,fileIndex}, Tfm,'bicubic','XData',[1 imgSize(2)],'YData',[1 imgSize(1)],'Size',imgSize));
    y   = I; 
    y = y / norm(y) ; 
    D(:,fileIndex) = y ;  
end

if para.saveStart
    save(fullfile(destDir, 'original-sralt.mat'),'D','xi_initial') ;
end

%% start the main loop
frOrig = cell(1,numImages) ;
T_in = cell(1,numImages) ;

T_ds = [ 0.5,   0, -0.5; ...
         0,   0.5, -0.5   ];
T_ds_hom = [ T_ds; [ 0 0 1 ]];

numIterOuter = 0 ; 
numIterInner = 0 ;

tic % time counting start
for scaleIndex = para.numScales:-1:1 
    
    iterNum = 0 ;  % iteration number of outer loop in each scale
    converged = 0 ;
    prevObj = inf ; % previous objective function value
    
    imgSize = para.canonicalImageSize / 2^(scaleIndex-1) ;    
    xi = cell(1,numImages) ;
    
    for fileIndex = 1 : numImages
             
        if scaleIndex == para.numScales
            T_in{fileIndex} = T_ds_hom^(scaleIndex-1)*transformations{fileIndex}*inv(T_ds_hom^(scaleIndex-1)) ;%T_ds_hom^(scaleIndex-1)单位矩阵 结果transformations{fileIndex}
        else
            T_in{fileIndex} = inv(T_ds_hom)*T_in{fileIndex}*T_ds_hom ;
        end
        
        % for display purposes
        if para.DISPLAY > 0
            fr = [1 1          imgSize(2) imgSize(2) 1; ...
                  1 imgSize(1) imgSize(1) 1          1; ...
                  1 1          1          1          1 ];
            % 
            frOrig{fileIndex} = T_in{fileIndex} * fr;
        end
        
    end
    
    while ~converged

        iterNum = iterNum + 1 ;
        numIterOuter = numIterOuter + 1 ;
        
        D= zeros(imgSize(1)*imgSize(2), numImages);
        J = cell(1,numImages) ;
        disp(['Scale ' num2str(scaleIndex) '  Iter ' num2str(iterNum)]) ;
        
        for fileIndex = 1 : numImages

            % transformed image and derivatives with respect to affine
            Tfm = fliptform(maketform('projective',T_in{fileIndex}'));

            I   = vec(imtransform(I0{scaleIndex,fileIndex}, Tfm,'bicubic','XData',[1 imgSize(2)],'YData',[1 imgSize(1)],'Size',imgSize));
            Iu  = vec(imtransform(I0x{scaleIndex,fileIndex},Tfm,'bicubic','XData',[1 imgSize(2)],'YData',[1 imgSize(1)],'Size',imgSize));
            Iv  = vec(imtransform(I0y{scaleIndex,fileIndex},Tfm,'bicubic','XData',[1 imgSize(2)],'YData',[1 imgSize(1)],'Size',imgSize));
            y   = I; %vec(I);

            Iu = (1/norm(y))*Iu - ( (y'*Iu)/(norm(y))^3 )*y ;
            Iv = (1/norm(y))*Iv - ( (y'*Iv)/(norm(y))^3 )*y ;

            y = y / norm(y) ; % normalize
            D(:,fileIndex) = y ;

            % transformation matrix to parameters;
            xi{fileIndex} = projective_matrix_to_parameters(para.transformType,T_in{fileIndex}) ;
            J{fileIndex} = image_Jaco(Iu, Iv, imgSize, para.transformType, xi{fileIndex});
        end
        
        lambda = 1*para.lambdac/sqrt(size(D,1)) ;
        
        % SRALT inner loop
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
        % using QR to orthogonalize the Jacobian matrix
        for fileIndex = 1 : numImages
            [Q{fileIndex}, R{fileIndex}] = qr(J{fileIndex},0) ;
        end
        
        D_tensor = refold(D',3,[imgSize(1) imgSize(2) numImages]);    

        if strcmp(para.optmet,'ADMM')
            [L, E, delta_xi, numIterInnerEach] = sralt_inner_admm(D_tensor, Q, lambda, para.alpha, para.p, para.inner_tol, para.inner_maxIter, para.muBound, para.rho0);
        end
        if strcmp(para.optmet,'ADMPG')
            [L, E, delta_xi, numIterInnerEach] = sralt_inner_admpg(D_tensor, Q, lambda, para.alpha, para.p, para.tau, para.inner_tol, para.inner_maxIter, para.muBound, para.rho0);
        end

        A = unfold(L,3)';
        E = unfold(E,3)';
        for fileIndex = 1 : numImages
            delta_xi{fileIndex} = inv(R{fileIndex})*delta_xi{fileIndex};
        end
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
        numIterInner = numIterInner + numIterInnerEach ;
        curObj = norm(svd(A),1) + lambda*norm(E(:),1) ;
        disp(['previous objective function: ' num2str(prevObj) ]);
        disp([' current objective function: ' num2str(curObj) ]);

        % step in paramters
        for i = 1 : numImages
            xi{i} = xi{i} + delta_xi{i};
            T_in{i} = parameters_to_projective_matrix(para.transformType,xi{i});
        end

        if para.DISPLAY > 0
            for i = 1 : numImages
                figure(1); clf ;
                imshow(I0{scaleIndex,i},[],'Border','tight');
                hold on;
                
                Tfm = fliptform(maketform('projective',inv(T_in{i}')));
                curFrame = tformfwd(fr(1:2,:)', Tfm )';
                plot( frOrig{i}(1,:),   frOrig{i}(2,:),   'g-', 'LineWidth', 2 );
                plot( curFrame(1,:), curFrame(2,:), 'r-', 'LineWidth', 2 );
                hold off;
                print('-f1', '-dbmp', fullfile(destDir, num2str(i))) ;
               
            end
        end
        
        if ( (prevObj - curObj < para.stoppingDelta) || iterNum >= para.maxIter )
            converged = 1;
            if ( prevObj - curObj >= para.stoppingDelta )
                disp('Maximum iterations reached') ;
            end
        else
            prevObj = curObj;
        end
        
    end
end

disp(['total number of iterations: ' num2str(numIterInner) ]);
disp(['number of outer loop: ' num2str(numIterOuter) ]);
timeConsumed = toc;

%% save the alignment results
Do = [] ; 
for fileIndex = 1 : numImages
    Tfm = fliptform(maketform('projective',T_in{fileIndex}'));   
    I   = vec(imtransform(I0{1,fileIndex}, Tfm,'bicubic','XData',[1 imgSize(2)],'YData',[1 imgSize(1)],'Size',imgSize));
    y   = I; 
    y = y / norm(y) ; % normalize
    Do = [Do y] ;
end

if para.saveEnd
    save(fullfile(destDir, 'final-sralt.mat'),'Do','A','E','xi') ;
    outputFileName = fullfile(destDir, 'final-sralt_results.txt'); 
    fid = fopen(outputFileName,'a') ;
    fprintf(fid, '%s\n', [' total number of iterations: ' num2str(numIterInner) ]) ;
    fprintf(fid, '%s\n', [' number of outer loop ' num2str(numIterOuter) ]) ;
    fprintf(fid, '%s\n', [' consumed time: ' num2str(timeConsumed)]) ;
    fprintf(fid, '%s\n', [' the parameters :']) ;
    fprintf(fid, '%s\n', [' transformType ' para.transformType ]) ;
    fprintf(fid, '%s\n', [' lambda ' num2str(para.lambdac) ' times sqrt(m)']) ;
    fprintf(fid, '%s\n', [' stoppingDelta of outer loop ' num2str(para.stoppingDelta) ]) ;
    fprintf(fid, '%s\n', [' stoppingDelta of inner loop ' num2str(para.inner_tol)]) ;
    if strcmp(para.optmet,'ADMM')
        fprintf(fid, '%s\n', [' optimization in inner loop is using ADMM algorithm']) ;
    else 
        fprintf(fid, '%s\n', [' optimization in inner loop is using LADM algorithm']) ;
    end
    fclose(fid);
end
close all;
