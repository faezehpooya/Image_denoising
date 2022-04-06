ORG_IM = imread('office_noisy.png'); % read the office noisy image
OFFICE_IM = imread('office.png'); % read the original image

A = im2double(ORG_IM);
AA = im2double(OFFICE_IM);

[m,n]=size(A); %size in pixels for A.
Anext=double(ORG_IM)/255; %use A=u(t) and Anext=u(t+h) (at the beginnig they are the same).

% print some initial objective quality evaluation MSE, SNR, and PSNR
err = immse(A,AA);
[peaksnr, snr] = psnr(A,AA);

fprintf('---------------------at First---------------------');
fprintf('\n The mean-squared error between office image and office noisy image is %0.4f\n', err);
fprintf('\n peak signal-to-noise ratio between office image and office noisy image is %0.4f\n',peaksnr);
fprintf(' \n The SNR value between office image and office noisy image is %0.4f\n', snr);

% setting the parameters

j= A;
j_plus_1 = zeros(size(j,1), size(j,2));


for iter = 1:250
    
    k= (1/iter)*0.25;
    north = zeros(size(j,1), size(j,2)); 
    north(2:end, 1:end) =  j(1:end-1, 1:end) ;
    north(1, :) = j(1, :); 
    
    dir_j_north = north - j;

    % South Gradient.
    south = zeros(size(j,1), size(j,2)); 
    south(1:end-1, 1:end) =  j(2:end, 1:end) ;
    south(end, :) = j(end, :); 

    dir_j_south = south - j;

    % West Gradient.
    west = zeros(size(j,1), size(j,2)); 
    west(:, 2:end) =  j(:, 1:end-1) ;
    west(:, 1) = j(:, 1); 

    dir_j_west = west - j;

    % East Gradient.
    east = zeros(size(j,1), size(j,2));
    east(:, 1:end-1) =  j(:, 2:end);
    east(:, end) = j(:, end); 

    dir_j_east = east - j;
    
    Dn = 1./(1+((dir_j_north./k).^2));
    Ds = 1./(1+((dir_j_south./k).^2));
    De = 1./(1+((dir_j_east./k).^2));
    Dw = 1./(1+((dir_j_west./k).^2));

    j_plus_1 = j + 0.25.*(Dn.*dir_j_north + Ds.*dir_j_south + De.*dir_j_east + Dw.*dir_j_west);
    j = j_plus_1; 
    
    %north = j(1:end-2, 2:end-1) - j(2:end-1, 2:end-1);
    %south = j(3:end, 2:end-1) - j(2:end-1, 2:end-1);
    %east = j(2:end-1, 3:end) - j(2:end-1, 2:end-1);
    %west = j(2:end-1, 1:end-2) - j(2:end-1, 2:end-1);

    %Dn = 1./(1+((north./k).^2));
    %Ds = 1./(1+((south./k).^2));
    %De = 1./(1+((east./k).^2));
    %Dw = 1./(1+((west./k).^2));

    %j_plus_1(2:end-1, 2:end-1) = j(2:end-1, 2:end-1) + 0.25.*(Dn.*north + Ds.*south + De.*east + Dw.*west);
    %j = j_plus_1; 
    
    
    if((iter==1)||(iter==5)||(iter==10)||(iter==30) ||(iter==100))
        figure; 
        subplot(1,2,1)
        imshow(im2double(ORG_IM));title('Office Noisy image');
        subplot(1,2,2)
        imshow(j);title(strcat('Office Noisy image after a diffusion time = ',int2str(iter)));
        pause(3);
        
  
        err = immse(j,AA);
        [peaksnr, snr] = psnr(j,AA);
        
        fprintf('---------------------in the for---------------------');
        fprintf('\n time is t = %d', iter);
        fprintf('\n The mean-squared error between office image and office noisy diffused image is %0.4f\n', err);
        fprintf('\n peak signal-to-noise ratio between office image and office noisy diffused image is %0.4f\n',peaksnr);
        fprintf('\n The SNR value between office image and office noisy diffused image is %0.4f\n', snr);
       
        KK = imabsdiff(im2double(ORG_IM),j);
        figure
        imshow(KK,[]);title(strcat('The imabsdiff difference between office noisy image and Smooth image at time t = ', int2str(iter)));
        pause(3);
    end
    
end

