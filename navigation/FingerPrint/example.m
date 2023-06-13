%% загрузка изображения
foldeName = "C:\Code\ipm-tasks\navigation\FingerPrint\img\";
imgName = "101_1.tif";
imgdata = imread(foldeName + imgName);

% отображение отпечатка
figure
imshow(imgdata);

%% обработка изображения

%% бинаризация изображения
clc
% использована простейшая бинаризация
% алгоритм Otsu - https://translated.turbopages.org/proxy_u/en-ru.ru.52ae5532-643d8278-e636a7a5-74722d776562/https/www.geeksforgeeks.org/binarization-of-digital-images-using-otsu-method-in-matlab/
%ad = binarize(imgdata);
%figure
%imshow(ad);
blocksize = 16;

bin = adaptiveThres(imgdata, blocksize);
figure
imshow(bin);

[p,z] = direction(bin,blocksize);
%% скелетизация
changing = true;
[rows, columns] = size(BW_Original);
BW_Thinned = BW_Original;
while changing
    BW_Del = ones(rows, columns);
    changing = false;
    % Setp 1
    for i=2:rows-1
        for j = 2:columns-1
            P = [BW_Thinned(i,j) BW_Thinned(i-1,j) BW_Thinned(i-1,j+1) BW_Thinned(i,j+1) BW_Thinned(i+1,j+1) BW_Thinned(i+1,j) BW_Thinned(i+1,j-1) BW_Thinned(i,j-1) BW_Thinned(i-1,j-1) BW_Thinned(i-1,j)]; % P1, P2, P3, ... , P8, P9, P2
            if (BW_Thinned(i,j) == 1 &&  sum(P(2:end-1))<=6 && sum(P(2:end-1)) >=2 && P(2)*P(4)*P(6)==0 && P(4)*P(6)*P(8)==0)   % conditions
                % No. of 0,1 patterns (transitions from 0 to 1) in the ordered sequence
                A = 0;
                for k = 2:size(P,2)-1
                    if P(k) == 0 && P(k+1)==1
                        A = A+1;
                    end
                end
                if (A==1)
                    BW_Del(i,j)=0;
                    changing = true;
                end
            end
        end
    end
    BW_Thinned = BW_Thinned.*BW_Del;  % the deletion must after all the pixels have been visited
    % Step 2 
    for i=2:rows-1
        for j = 2:columns-1
            P = [BW_Thinned(i,j) BW_Thinned(i-1,j) BW_Thinned(i-1,j+1) BW_Thinned(i,j+1) BW_Thinned(i+1,j+1) BW_Thinned(i+1,j) BW_Thinned(i+1,j-1) BW_Thinned(i,j-1) BW_Thinned(i-1,j-1) BW_Thinned(i-1,j)];
            if (BW_Thinned(i,j) == 1 && sum(P(2:end-1))<=6 && sum(P(2:end-1)) >=2 && P(2)*P(4)*P(8)==0 && P(2)*P(6)*P(8)==0)   % conditions
                A = 0;
                for k = 2:size(P,2)-1
                    if P(k) == 0 && P(k+1)==1
                        A = A+1;
                    end
                end
                if (A==1)
                    BW_Del(i,j)=0;
                    changing = true;
                end
            end
        end
    end
    BW_Thinned = BW_Thinned.*BW_Del;
end%while

figure
subplot(1,2,1)
imshow(BW_Original)
subplot(1,2,2)
imshow(BW_Thinned)

%% выделение особых точек
BW_morph = tmpDelete(BW_Original); 
figure
subplot(1,2,1)
imshow(BW_Original)
subplot(1,2,2)
imshow(BW_morph)
%% функции
function binImg = binarize(img)
    binImg = [];
    for i = 1 : size(img,1)
        tmp = [];
        for j = 1 : size(img,2 )
            pixel = img(i,j);
            if pixel > 100
                pixel = 1;
            else
                pixel = 0;
            end
            binImg(i, j) = pixel;
        end
    end
end

function [final] = fftenhance(image,f)
    I = 255-double(image);
    [w,h] = size(I);
    %out = I;
    w1=floor(w/32)*32;
    h1=floor(h/32)*32;
    inner = zeros(w1,h1);
    for i=1:32:w1
        for j=1:32:h1
            a=i+31;
            b=j+31;
            F=fft2( I(i:a,j:b) );
            factor=abs(F).^f;
            block = abs(ifft2(F.*factor));
            larv=max(block(:));
            if larv==0
                larv=1;
            end
            block= block./larv;
            inner(i:a,j:b) = block;
        end
    end
    final=uint8(inner*255);
    %final=histeq(final);
end

function [o] = adaptiveThres(a,W)
    %Adaptive thresholding is performed by segmenting image a
    [w,h] = size(a);
    o = zeros(w,h);
    %seperate it to W block
    %step to w with step length W
    for i=1:W:w
        for j=1:W:h
            mean_thres = 0;
            if i+W-1 <= w && j+W-1 <= h
                mean_thres = mean2(a(i:i+W-1, j:j+W-1));
                mean_thres = 0.8*mean_thres;
                o(i:i+W-1,j:j+W-1) = a(i:i+W-1,j:j+W-1) < mean_thres;
            end
        end
    end
end

function [p,z] = direction(image,blocksize)
    %image=adaptiveThres(image,16,0);
    [w,h] = size(image);
    direct = zeros(w,h);
    gradient_times_value = zeros(w,h);
    gradient_sq_minus_value = zeros(w,h);
    gradient_for_bg_under = zeros(w,h);
    W = blocksize;
    theta = 0;
    sum_value = 1;
    bg_certainty = 0;
    blockIndex = zeros(ceil(w/W),ceil(h/W));
    %directionIndex = zeros(ceil(w/W),ceil(h/W));
    times_value = 0;
    minus_value = 0;
    center = [];
    filter_gradient = fspecial('sobel');
    %to get x gradient
    I_horizontal = filter2(filter_gradient,image);
    %to get y gradient
    filter_gradient = transpose(filter_gradient);
    I_vertical = filter2(filter_gradient,image);
    gradient_times_value=I_horizontal.*I_vertical;
    gradient_sq_minus_value=(I_verticalI_horizontal).*(I_vertical+I_horizontal);
    gradient_for_bg_under = (I_horizontal.*I_horizontal) + (I_vertical.*I_vertical);
    for i=1:W:w
        for j=1:W:h
            if j+W-1 < h && i+W-1 < w
                times_value = sum(sum(gradient_times_value(i:i+W-1, j:j+W-1)));
                minus_value = sum(sum(gradient_sq_minus_value(i:i+W-1, j:j+W1)));
                sum_value = sum(sum(gradient_for_bg_under(i:i+W-1, j:j+W-1)));
                bg_certainty = 0;
                theta = 0;
    
                if sum_value ~= 0 && times_value ~=0
                    %if sum_value ~= 0 & minus_value ~= 0 & times_value ~= 0
                    bg_certainty = (times_value*times_value + minus_value*minus_value)/(W*W*sum_value);
                    if bg_certainty > 0.05
                        blockIndex(ceil(i/W),ceil(j/W)) = 1;
                        %tan_value = atan2(minus_value,2*times_value);
                        tan_value = atan2(2*times_value,minus_value);
                        theta = (tan_value)/2 ;
                        theta = theta+pi/2;
                        center = [center;[round(i + (W-1)/2),round(j + (W-1)/2),theta]];
                    end
                end
            end
            times_value = 0;
            minus_value = 0;
            sum_value = 0;
        end
    end
    imagesc(direct);
    hold on;
    [u,v] = pol2cart(center(:,3),8);
    quiver(center(:,2),center(:,1),u,v,0,'g');
    hold off;
    x = bwlabel(blockIndex,4);
    y = bwmorph(x,'close');
    z = bwmorph(y,'open');
    p = bwperim(z);
end