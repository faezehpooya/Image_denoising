ORG_IM = imread('office_noisy.png'); % read the office noisy image
Sigma_values = {0.5,1,2,5,10,50};

for j = 1:length(Sigma_values)% for each of them sigma values 
    sigma = Sigma_values{j};
    h = fspecial('gaussian',[3*ceil(sigma) 3*ceil(sigma)],sigma); % create the [3*sigma 3*sigma] Gaussian kernel  
    %Im_smooth = imfilter(ORG_IM,h); % perform the the kernel on the image  
    I_smooth = imgaussfilt(im2double(ORG_IM),sigma);
    
    % plotting two images
    figure;   
    subplot(1,2,1);imshow(im2double(ORG_IM),[]);title('Noisy image');
    subplot(1,2,2);imshow(I_smooth,[]);title(strcat('Smooth image with sigma = ', num2str(sigma)));
    pause(3);
    
    % plotting the difference between images
    K = imabsdiff(im2double(ORG_IM),I_smooth);
    figure
    imshow(K,[]);title(strcat('The imabsdiff difference between office noisy image and Smooth image with sigma = ', num2str(sigma)));
    pause(3);
   
end


