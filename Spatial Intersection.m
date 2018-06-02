% IMBDC Exercise 4 Spatial intersection
% Worked by Karam Mawas 2946939

close all
clear all
clc

% points are already measured by the ginput function
% x = [135.437062937063,3565.83706293706,3955.65524475524,3907.67762237762,3283.96853146853,2270.44125874126,2612.28181818182,1850.63706293706,2426.36853146853]';
% y = [177.081818181819,2348.06923076923,2072.19790209790,1556.43846153846,650.860839160840,1676.38251748252,1220.59510489511,428.964335664336,380.986713286714]';
% Note
% R0020774=im1\R0020813=im2\R0020814=im3\R0020815=im4\R0020816=im5\R0020849=im6\R0020850=im7\R0020851=im8\R0020852=im9

% reading the images
% im71 = imread('C:\Users\ehteshamhasnain\Downloads\Karam\Exercise 2\Ricoh-UAV_Subblock\R0020771.jpg');
% im72 = imread('C:\Users\ehteshamhasnain\Downloads\Karam\Exercise 2\Ricoh-UAV_Subblock\R0020772.jpg');
% im73 = imread('C:\Users\ehteshamhasnain\Downloads\Karam\Exercise 2\Ricoh-UAV_Subblock\R0020773.jpg');
% im74 = imread('C:\Users\ehteshamhasnain\Downloads\Karam\Exercise 2\Ricoh-UAV_Subblock\R0020774.jpg');
% im75 = imread('C:\Users\ehteshamhasnain\Downloads\Karam\Exercise 2\Ricoh-UAV_Subblock\R0020775.jpg');
% im12 = imread('C:\Users\ehteshamhasnain\Downloads\Karam\Exercise 2\Ricoh-UAV_Subblock\R0020812.jpg');
% im13 = imread('C:\Users\ehteshamhasnain\Downloads\Karam\Exercise 2\Ricoh-UAV_Subblock\R0020813.jpg');
% im14 = imread('C:\Users\ehteshamhasnain\Downloads\Karam\Exercise 2\Ricoh-UAV_Subblock\R0020814.jpg');
% im15 = imread('C:\Users\ehteshamhasnain\Downloads\Karam\Exercise 2\Ricoh-UAV_Subblock\R0020815.jpg');
% im16 = imread('C:\Users\ehteshamhasnain\Downloads\Karam\Exercise 2\Ricoh-UAV_Subblock\R0020816.jpg');
% im49 = imread('C:\Users\ehteshamhasnain\Downloads\Karam\Exercise 2\Ricoh-UAV_Subblock\R0020849.jpg');
% im50 = imread('C:\Users\ehteshamhasnain\Downloads\Karam\Exercise 2\Ricoh-UAV_Subblock\R0020850.jpg');
% im51 = imread('C:\Users\ehteshamhasnain\Downloads\Karam\Exercise 2\Ricoh-UAV_Subblock\R0020851.jpg');
% im52 = imread('C:\Users\ehteshamhasnain\Downloads\Karam\Exercise 2\Ricoh-UAV_Subblock\R0020852.jpg');
% im53 = imread('C:\Users\ehteshamhasnain\Downloads\Karam\Exercise 2\Ricoh-UAV_Subblock\R0020853.jpg');
% im71 = rgb2gray(im71); im72 = rgb2gray(im72);

% initializing the interior orientatin of the camera and rotation matricies
% for all the images

% Interior Orientation for the images (the K matrix) 
K  = [-3933.36363636,-0.00000000,2143.50000000;
      -0.00000000,3933.36363636,1423.50000000;
      0.00000000,0.00000000,1.00000000];
% rotation matricies

% R_71 = [-0.888776619260	 -0.457304048725	 -0.030807922295;
% 	  0.458330362120	 -0.886297121687	 -0.066413035232;
% 	  0.003065977044	 -0.073146559110	  0.997316489724];
% X0_71 =  [512906.565340;	  5427602.989520;	      515.744690];
% 
% R_72 = [-0.870960174003	 -0.488580136412	 -0.052132769008;
% 	  0.488253956273	 -0.872479485577	  0.019688103784;
% 	 -0.055103987919	 -0.008306476422	  0.998446069132];
% X0_72 =  [512891.972680;	  5427626.881650;	      515.419150];
%   
% R_73 = [-0.853199280893	 -0.519667012797	  0.044689852260;
% 	  0.517892492096	 -0.854220979117	 -0.045758993294;
% 	  0.061954448707	 -0.015897001214	  0.997952369423];
% X0_73 =  [512876.458510;	  5427651.545480;	      514.346390];
  
R_74 = [-0.898364354805	 -0.424038535386	 -0.114598457768;
	  0.431222964600	 -0.901052927842	 -0.046372147116;
	 -0.083595698552	 -0.091076570722	  0.992328986500];  
X0_74= [512862.179910;	  5427675.837710;	      514.907080];  
  
% R_75 = [ -0.865048359102	 -0.501663541135	  0.005002790291;
% 	  0.500877525330	 -0.864173861893	 -0.048220753218;
% 	  0.028513874422	 -0.039207498226	  0.998824174241];
% X0_75 = [512846.973530;	  5427700.537190;	      515.308980];
%   
% R_12 = [ 0.886120155244	  0.448624945312	 -0.116304466445;
% 	 -0.443803228369	  0.893694064736	  0.065951596980;
% 	  0.133528142953	 -0.006824741672	  0.991021522441];  
% X0_12= [512892.209370;	  5427714.738940;	      515.158090];    
  
R_13 = [0.872990632012,	  0.487238787937,	  0.022039055070;
    -0.486140187346,	  0.872897199467,	 -0.041451144856;
    -0.039434435028,	  0.025472390784,	  0.998897433494];
X0_13 =  [512906.392000;	  5427691.184110;	      516.570080];

R_14 = [0.867281138830	  0.493392418024	 -0.066237059618;
	 -0.493411639320	  0.869625609377	  0.017212021948;
	  0.066093724462	  0.017754474173	  0.997655450661];
X0_14 =  [512920.853830;	  5427668.335460;	      516.113920];
  
R_15 = [0.858414149205	  0.508703759544	 -0.065921418897;
	 -0.507788494956	  0.860915774482	  0.031222966529;
	  0.072636029862	  0.006671901839	  0.997336198527];
X0_15 =  [512936.566090;	  5427645.080020;	      515.434980];
  
R_16 = [0.854735681670	  0.514751542165	  0.066766490992;
	 -0.512025217779	  0.857252179862	 -0.054303558625;
	 -0.085188560491	  0.012229062108	  0.996289796797];  
X0_16 = [512952.422120;	  5427622.063810;	      515.982850];
  
R_49 = [-0.867233884963	 -0.496431897214	 -0.038219892721;
	  0.497756374147	 -0.862573995849	 -0.090579764190;
	  0.011999198600	 -0.097578036019	  0.995155538657];  
X0_49 = [513010.498970;	  5427654.064230;	      514.392000];

R_50 = [-0.839539757811	 -0.528716554759	 -0.125027196158;
	  0.538291662815	 -0.840647830782	 -0.059609649767;
	 -0.073587232584	 -0.117345768246	  0.990360989678];  
X0_50 = [512996.502990;	  5427678.079250;	      513.677870];
  
R_51 = [-0.827570749773	 -0.557530411593	  0.065471324026;
	  0.549857636160	 -0.828569770064	 -0.105492730048;
	  0.113062965097	 -0.051302790236	  0.992262460056];  
X0_51 = [512980.99511;	  5427701.52671;	      514.79429];

R_52 = [-0.866544118547	 -0.490057510421	 -0.094577624687;
	  0.490759244246	 -0.871123481449	  0.017298677877;
	 -0.090866136699	 -0.031424776041	  0.995367182829];  
X0_52 = [512966.499230;	  5427725.063170;	      517.490680];
  
% R_53 = [	-0.857743949672	 -0.514007847249	 -0.008440957694;
% 	  0.513687553089	 -0.857617486843	  0.024846409457;
% 	 -0.020010362361	  0.016975842479	  0.999655643795];  
% X0_13 = [512951.374930;	  5427750.327570;	      515.903230];

% to be consistent with my type of naming

R1 = R_74; R2 = R_13; R3 = R_14; R4 = R_15; R5 = R_16; R6 = R_49;
R7 = R_50; R8 = R_51; R9 = R_52;
X0_1 = X0_74; X0_2 = X0_13; X0_3 = X0_14; X0_4 = X0_15; X0_5 = X0_16;
X0_6 = X0_49; X0_7 = X0_50; X0_8 = X0_51; X0_9 = X0_52;

% Computing the Projection matricies
P1 = K*R1*[eye(3) -X0_1]; 
P2 = K*R2*[eye(3) -X0_2];
P3 = K*R3*[eye(3) -X0_3];
P4 = K*R4*[eye(3) -X0_4];
P5 = K*R5*[eye(3) -X0_5];
P6 = K*R6*[eye(3) -X0_6];
P7 = K*R7*[eye(3) -X0_7]; 
P8 = K*R8*[eye(3) -X0_8]; 
P9 = K*R9*[eye(3) -X0_9];


x_h = ones(2,9); 
x = ones(1,9); y = ones(1,9); f = ones(1,9);
for i=1:9
    m = imread(['E:\UNI\UNI\Third semester\Image-based Data Collection\Labs\Exercise 4\im' num2str(i) '.jpg']);
    f(i) = figure(i);
%    imshow(rgb2gray(m))
    imshow(m)
    [x(i),y(i)] = ginput(1);
    % x = imtool(rgb2gray(m));
    % [x,y] = getline('close')
    hold on; plot(x(i),y(i),'r+');

    % The homogeneous coordinates of the images points
    x_h(1,i) = x(i); x_h(2,i) = y(i); 
        
end

% computing the A matrix
A = [
    x_h(1,1)*P1(3,:) - P1(1,:)
    x_h(2,1)*P1(3,:) - P1(2,:)
    x_h(1,2)*P2(3,:) - P2(1,:)
    x_h(2,2)*P2(3,:) - P2(2,:)
    x_h(1,3)*P3(3,:) - P3(1,:)
    x_h(2,3)*P3(3,:) - P3(2,:)
    x_h(1,4)*P4(3,:) - P4(1,:)
    x_h(2,4)*P4(3,:) - P4(2,:)
    x_h(1,5)*P5(3,:) - P5(1,:)
    x_h(2,5)*P5(3,:) - P5(2,:)
    x_h(1,6)*P6(3,:) - P6(1,:)
    x_h(2,6)*P6(3,:) - P6(2,:)
    x_h(1,7)*P7(3,:) - P7(1,:)
    x_h(2,7)*P7(3,:) - P7(2,:)
    x_h(1,8)*P8(3,:) - P8(1,:)
    x_h(2,8)*P8(3,:) - P8(2,:)
    x_h(1,9)*P9(3,:) - P9(1,:)
    x_h(2,9)*P9(3,:) - P9(2,:)
     ];
 
 % Extracting the Object coordinates of  the points
 [U, V, X] = svd(A,0);
 X = X(:,end);
 
 % Coordinates of 3d point after normalisation
 Terrain_Point_3D = X(:)./X(4);
 
 
 % computing the back transformation for the points 
 % the normalized images coordinates 
 x_trafo_1 = P1*Terrain_Point_3D; x_trafo_1 = x_trafo_1(:)./x_trafo_1(3);
 x_trafo_2 = P2*Terrain_Point_3D; x_trafo_2 = x_trafo_2(:)./x_trafo_2(3);
 x_trafo_3 = P3*Terrain_Point_3D; x_trafo_3 = x_trafo_3(:)./x_trafo_3(3);
 x_trafo_4 = P4*Terrain_Point_3D; x_trafo_4 = x_trafo_4(:)./x_trafo_4(3);
 x_trafo_5 = P5*Terrain_Point_3D; x_trafo_5 = x_trafo_5(:)./x_trafo_5(3);
 x_trafo_6 = P6*Terrain_Point_3D; x_trafo_6 = x_trafo_6(:)./x_trafo_6(3);
 x_trafo_7 = P7*Terrain_Point_3D; x_trafo_7 = x_trafo_7(:)./x_trafo_7(3);
 x_trafo_8 = P8*Terrain_Point_3D; x_trafo_8 = x_trafo_8(:)./x_trafo_8(3);
 x_trafo_9 = P9*Terrain_Point_3D; x_trafo_9 = x_trafo_9(:)./x_trafo_9(3);
 x_trafo = [x_trafo_1,x_trafo_2,x_trafo_3,x_trafo_4,x_trafo_5,x_trafo_6,x_trafo_7,x_trafo_8,x_trafo_9];
 
 % plotting the computed points
 for i=1:9
     m = imread(['E:\UNI\UNI\Third semester\Image-based Data Collection\Labs\Exercise 4\im' num2str(i) '.jpg']);
    f(i) = figure(i);
    imshow(m)
    plot(x_h(1,i), x_h(2,i), 'r+'); hold on; 
    plot(x_trafo(1,i), x_trafo(2,i),'b+')
 end
 
 % difference between the computed image coordinates and the original ones
 diff = ones(2,9); diff_x = ones(1,9); diff_y = ones(1,9); 
 for i=1:9
     for j = 1:2
         diff(j,i) = x_h(j,i) - x_trafo(j,i);
     end
     diff_x(1,i) = diff(1,i)^2;
     diff_y(1,i) = diff(2,i)^2;
 end
 
 % Error Estimation
 VV_x = sum(diff_x,2); VV_y = sum(diff_y,2);
 sigma_0_x = sqrt(VV_x/(2*9-3))
 sigma_0_y = sqrt(VV_y/(2*9-3))
 
 