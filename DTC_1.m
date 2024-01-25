%% Diurnal Temperature Cycle Models
%To estimate LST with hourly sequence, a model consisting of a harmonic and an exponential term is fitted to the Earths DTC model, which describe the effect of the sun and the reduction of LST at night, respectively. Modeling estimates the parameters that describe DTC and can also be useful for interpolating missing data due to technical or cloud problems. The parameters depend on all modeled temperatures and are therefore hardly affected by outliers.
%DTC:: A two-part, semi-empirical DTC model is developed in [Gottsche, 2001]. The LST cycle is modeled using the cosine function to predict the evolution of LST during the day and based on the thermal diffusion equation and an exponential function to describe the temperature decrease at night, assuming that natural surfaces follow Newton s law of cooling

%!! Important Note!! 
%In this code, 4 MODIS products are used for the daily temperature cycle model.
%which measure the land surface temperature at 10:30 AM/PM local solar time for the Terra satellite and 1:30 AM/PM for the Aqua satellite and are freely available as MOD11A1 and MYD11A1.
%To run this code, it is necessary to store the LST images of each of these
%four hours in a separate folder and number them in some way so that the image of each date can be read at the same time, for example, if we want to do modeling for a year. We have 365 images in each folder, which are numbered from 001 to 365.
% This code can estimate the land surface temperature image in 24 hours by having four MODIS LST images.

%%!!Code execution steps
% 1-Change LST folder paths in line 27, 28, 29, 30, 150, 160.
% 2-Enter the Latitude and longitude of your study area in line 83, 100.
% 3-Enter the number i based on the number of input images in line 33.
% 4- In the line 117, enter the product of the number of rows and columns of your image against i.
% 5-Enter the dimensions of your image in line 147, 157.

%% This code was written by Fahime Arabi. if you have any questions about it, I will answer you with the following email:
%Fahimearabi1993@gmail.com

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
clc;
clear;
close all;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
path_LST1='F:\Thesis\5-LST\modis\lst\2020\1_LST_SSA_MOD\2-Night'
path_LST10='F:\Thesis\5-LST\modis\lst\2020\1_LST_SSA_MOD\1-Day'
path_LST13='F:\Thesis\5-LST\modis\lst\2020\2-LST_SSA_MYD\3-DAY_TIFF'
path_LST22='F:\Thesis\5-LST\modis\lst\2020\2-LST_SSA_MYD\4-NIGHT_TIFF'
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
 %Enter the number i based on the number of input images.
for i=1:365
    
    cd(path_LST1);
    a=dir('*.tif');
    nameLST1=a(i).name;
    [img,Ref]=geotiffread([path_LST1, '\', nameLST1]);
    strname=strsplit(nameLST1,'_');
    ddatte=strname{1,2};
    info=geotiffinfo([path_LST1, '\', nameLST1]);
    img=double(img);
    LST1=(img(:))';
    LST1=LST1*0.02;
% %     0.02 is a coefficient that must be multiplied in MODIS LST images
% to convert to Kelvin.

    cd(path_LST10);
    a=dir('*.tif');
    nameLST10=a(i).name;
    [img,Ref]=geotiffread([path_LST10, '\', nameLST10]);
     [B4,RefMatrx]=geotiffread(nameLST10);
    InfoB4=geotiffinfo([path_LST10, '\', nameLST10]);
    img=double(img);
    LST10=(img(:))';
    LST10=LST10*0.02;
    
    
    cd(path_LST13);
    a=dir('*.tif');
    nameLST13=a(i).name;
    [img,Ref]=geotiffread([path_LST13, '\', nameLST13]);
    info=geotiffinfo([path_LST13, '\', nameLST13]);
    img=double(img);
    LST13=(img(:))';
    LST13=LST13*0.02;
    
    
    cd(path_LST22);
    a=dir('*.tif');
    nameLST22=a(i).name;
    [img,Ref]=geotiffread([path_LST22, '\', nameLST22]);
    info=geotiffinfo([path_LST22, '\', nameLST22]);
    img=double(img);
    LST22=(img(:))';
    LST22=LST22*0.02;
    
 data=cat(1,LST10,LST13,LST22,LST1);%Image Read
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%In this section, sunrise, sunset and day length are calculated.
    d=i; 
    %Enter the Latitude of your study area
    Latitude=31.0833;
    declination=23.45*sin((2*pi*(d-80))/(365));
 
    Axis=23.439*pi/180;
    j=pi/182.625;
    m=1-tan(Latitude*pi/180).*tan(Axis*cos(j*d));

     m(m>2)=2;
     m(m<0)=0;
     b=acos(1-m)/pi;
     hours=b*24;
     w=hours;  %Length of day
     sunrise=12-(w/2); %Sunrise
     sunset=sunrise+w;  %Sunset
     ts=sunset- sunrise;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%In this section, local noon time is calculated
 longitude=54.3569;
%%Enter the longitude of your study area
    delta_UTC=3.5;
    d=i;
    LT=12;
    LSTM=15*delta_UTC;
    DC=2*pi/365;
    B=DC*(d+10)+0.033*sin(DC*(d-2));
    E0T=(9.87*(sin(2*B))+(7.6)*(sin(B-0.2)));
    TC=4*(longitude-LSTM)+E0T;
    LST=LT+(TC/60);
    tm=LST+1; %tm
    tm=tm- sunrise;
    u=(pi/w)*(ts-tm);
    k=(w/pi)*atan((pi/w)*(ts-tm))
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the line below, enter the product of the number of rows and columns of your image against i.
 for i= 1:115116 
    y=[data(:,i)];
    f=[1, cos(((pi/w)*((10- sunrise)-tm))); 1, cos(((pi/w)*((13- sunrise)-tm)));1, cos(((pi/w)*((22- sunrise)-tm)))*exp(-((22- sunrise)-ts)/k);1, cos(((pi/w)*((25- sunrise)-tm)))*exp(-((25- sunrise)-ts)/k)];;               
    a=(inv(transpose(f)*f))*transpose(f)*double(y);
    result(:,i)=a;
end
mean=result(1,:);
amplitude=result(2,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tmax=mean+amplitude*cos((pi/w)*(tm-tm));
Tsrd=mean+amplitude*cos((pi/w)*((sunrise-1)-(LST+1)));
Tsrd1 =(mean)+((amplitude*cos((pi/w)*(sunset-(LST+1))))*exp(-(24-sunset)/k));

Ta=(Tmax-Tsrd)/(cos(pi/4)+1);
T0=Tsrd+Ta.*cos(pi/4);
omrgaT=Ta.*(cos(u).*(Ta*cos(u)+T0-Tsrd1)+((pi/4).*sin(u).*(24-ts)*(T0-Tsrd1)))/(T0-Tsrd1-Ta.*((pi/4).*sin(u).*(24-ts)-cos(u)));

u=(pi/w)*(ts-tm);


 u=(pi/w)*(ts-tm);

k=((Ta.*(cos(u))-omrgaT))/((Ta.*pi/w)*sin(u));
ts=round(ts)
tss=ts+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:ts
    t=i;
    T =mean+amplitude.*cos((pi/w).*((t)-tm));
  % Enter the dimensions of your image
    T=reshape(T,318,362); 
   nnaammee=['lst_' num2str(ddatte) '_' num2str(i)];
  % Change LST folder paths.
   filenam = ['F:\Thesis\4-LST\6-DCY\YEARS\hourly_size\' nnaammee '.tif']
   geotiffwrite( filenam,T,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);
end
for i=tss:24
    t=i;
    T =(mean+omrgaT)+((amplitude*cos((pi/w)*(ts-tm))-omrgaT)*exp(-(t-ts)/k));
     % Enter the dimensions of your image
    T=reshape(T,318,362); 
    nnaammee=['lst_' num2str(ddatte) '_' num2str(i)];
    %Change LST folder paths.
   filenam = ['F:\Thesis\4-LST\6-DCY\YEARS\hourly_size\' nnaammee '.tif']
   geotiffwrite( filenam,T,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);
 
end
 end
