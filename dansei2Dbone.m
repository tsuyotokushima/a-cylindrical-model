% clear
% close all
% clc
% 
% %%
% velocity_longitudinal_water=1500; %水中縦波音速
% velocity_share_water=0; %水中横波音速
% frequency=1.0e+6;
% rhow=1000; %水密度
% dimension=2; %次元
% Kwl=rhow*velocity_longitudinal_water^2; %縦波　定数
% 
% rhob=2000; %骨密度
% velocity_longitudinal_bone=3500; %骨中縦波音速
% velocity_share_bone=2400; %骨中横波音速
% 
% dx=velocity_longitudinal_water/frequency/40; %空間分解能
% dy=dx;
% dt=dx/velocity_longitudinal_bone/sqrt(dimension); %時間分解能
% 
% model=zeros(500,500);%シミュレーション空間
% nt=600; %for文回す回数
% 
% nx=300;
% ny=100;
% 
% %%
% %入力波形
% n=round(1/frequency/dt);
% wave(:)=sin(2*pi*frequency*dt*(1:n));
% % wd(:)=0.5-0.5*cos(2*pi*frequency*dt*(1:n));
% % wave2=wave.*wd;
% % 
% % figure;plot(wave2)
% 
% %%
% Txx(:,:)=single(zeros(nx,ny));
% Tyy(:,:)=single(zeros(nx,ny));
% Txy(:,:)=single(zeros(nx+1,ny+1));
% Ux(:,:)=single(zeros(nx+1,ny));
% Uy(:,:)=single(zeros(nx,ny+1));
% 
% C11(:,:)=single(zeros(nx,ny)); %Txx依存
% C22(:,:)=single(zeros(nx,ny)); %Tyy依存
% C66(:,:)=single(zeros(nx+1,ny+1)); %Txy依存
% C12(:,:)=single(zeros(nx,ny)); %Txx依存
% 
% 
% lowx(:,:)=single(zeros(nx+1,ny)); %Ux依存　たてなみ関連
% lowy(:,:)=single(zeros(nx,ny+1)); %Uy依存　よこなみ関連
% 
% %%
% %全て水
% lowx(:,:)=rhow;
% lowy(:,:)=rhow;
% C11(:,:)=rhow*velocity_longitudinal_water^2;
% C22(:,:)=rhow*velocity_longitudinal_water^2;
% C66(:,:)=rhow*velocity_share_water^2;
% 
% C12(:,:)=rhow*velocity_longitudinal_water^2-2*rhow*velocity_share_water^2;
% C21(:,:)=C12(:,:);
% %C12(:,:)=C11-2*C66;
% 
% 
% %%
% %骨ゾーン   ブロック形　　上下左右に水の層
% lowx(102:200,20:80)=rhob;
% lowy(101:200,22:80)=rhob;
% C11(101:200,20:80)=rhob*velocity_longitudinal_bone^2;
% C22(101:200,20:80)=rhob*velocity_longitudinal_bone^2;
% C66(102:200,22:80)=rhob*velocity_share_bone^2;
% 
% C12(101:200,20:80)=rhob*velocity_longitudinal_bone^2-2*rhob*velocity_share_bone^2;
% %C12=C11-2*C66;
% C21(:,:)=C12(:,:);
% 
% %%
% %中間ゾーン
% lowx(101,20:80)=(rhow+rhob)/2;
% lowx(201,20:80)=(rhow+rhob)/2;
% lowy(101:200,21)=(rhow+rhob)/2;
% lowy(101:200,81)=(rhow+rhob)/2;
% 
% inputpointx=1;%音波入力位置　
% inputpointy=single(ny/2-48:ny/2+48);%音波入力位置
% 
% figure;
% %% FDTD
% for s=1:nt
%     if s<=n
%         Txx(inputpointx,inputpointy)=wave(s);
%         Tyy(inputpointx,inputpointy)=wave(s);
%     else
%         Txx(inputpointx,inputpointy)=0;
%         Tyy(inputpointx,inputpointy)=0;
%     end
%     
% %FDTDの式 
%         
%         Ux(2:nx,1:ny)=Ux(2:nx,1:ny)+dt./lowx(2:nx,1:ny).*((Txx(2:nx,1:ny)-Txx(1:nx-1,1:ny))/dx+(Txy(2:nx,2:ny+1)-Txy(2:nx,1:ny))/dy);
%         Uy(1:nx,2:ny)=Uy(1:nx,2:ny)+dt./lowy(1:nx,2:ny).*((Tyy(1:nx,2:ny)-Tyy(1:nx,1:ny-1))/dy+(Txy(2:nx+1,2:ny)-Txy(1:nx,2:ny))/dx);
%     
%         Txx(1:nx,1:ny)=Txx(1:nx,1:ny)+dt.*(C11(1:nx,1:ny).*(Ux(2:nx+1,1:ny)-Ux(1:nx,1:ny))/dx+C12(1:nx,1:ny).*(Uy(1:nx,2:ny+1)-Uy(1:nx,1:ny))/dy);
%         Tyy(1:nx,1:ny)=Tyy(1:nx,1:ny)+dt.*(C21(1:nx,1:ny).*(Ux(2:nx+1,1:ny)-Ux(1:nx,1:ny))/dx+C22(1:nx,1:ny).*(Uy(1:nx,2:ny+1)-Uy(1:nx,1:ny))/dy);
%         Txy(2:nx,2:ny)=Txy(2:nx,2:ny)+dt.*C66(2:nx,2:ny).*((Ux(2:nx,2:ny)-Ux(2:nx,1:ny-1))/dy+(Uy(2:nx,2:ny)-Uy(1:nx-1,2:ny))/dx);
%         
% 
%        
%     
%     T=sqrt(Txx.^2+Tyy.^2);
%     imagesc(T)
%     axis equal
% 
%     M=getframe
%     
%     
%     clc
%     s
% end



clear
close all
clc

%%
velocity_longitudinal_water=1500; %水中縦波音速
velocity_share_water=0; %水中横波音速
frequency=1.0e+6;
rhow=1000; %水密度
dimension=2; %次元
Kwl=rhow*velocity_longitudinal_water^2; %縦波　定数

rhob=2000; %骨密度
velocity_longitudinal_bone=3500; %骨中縦波音速
velocity_share_bone=2400; %骨中横波音速

dx=velocity_longitudinal_water/frequency/40; %空間分解能
dy=dx;
dt=dx/velocity_longitudinal_water/sqrt(dimension); %時間分解能

model=zeros(500,500);%シミュレーション空間
nt=100; %for文回す回数

nx=300;
ny=100;

%%
%入力波形
n=round(1/frequency/dt);
wave(:)=sin(2*pi*frequency*dt*(1:n));
% wd(:)=0.5-0.5*cos(2*pi*frequency*dt*(1:n));
% wave2=wave.*wd;
% 
% figure;plot(wave2)

%% ALLOCATE MEMORY
Txx(:,:)=single(zeros(nx,ny));
Tyy(:,:)=single(zeros(nx,ny));
Txy(:,:)=single(zeros(nx+1,ny+1));
Ux(:,:)=single(zeros(nx+1,ny+1));
Uy(:,:)=single(zeros(nx+1,ny+1));

C11(:,:)=single(zeros(nx,ny)); %Txx依存
C22(:,:)=single(zeros(nx,ny)); %Tyy依存
C66(:,:)=single(zeros(nx+1,ny+1)); %Txy依存
C12(:,:)=single(zeros(nx,ny)); %Txx依存


lowx(:,:)=single(zeros(nx+1,ny)); %Ux依存　たてなみ関連
lowy(:,:)=single(zeros(nx,ny+1)); %Uy依存　よこなみ関連

%% Coefficients for water
%全て水
lowx(:,:)=rhow;
lowy(:,:)=rhow;
C11(:,:)=rhow*velocity_longitudinal_water^2;
C22(:,:)=rhow*velocity_longitudinal_water^2;
C66(:,:)=rhow*velocity_share_water^2;

C12(:,:)=rhow*velocity_longitudinal_water^2-2*rhow*velocity_share_water^2;
C21(:,:)=C12(:,:);
%C12(:,:)=C11-2*C66;


%% Coefficients for bone
% %骨ゾーン   ブロック形　　上下左右に水の層
% lowx(52:200,20:80)=rhob;
% lowy(51:200,22:80)=rhob;
% C11(51:200,20:80)=rhob*velocity_longitudinal_bone^2;
% C22(51:200,20:80)=rhob*velocity_longitudinal_bone^2;
% C66(52:200,22:80)=rhob*velocity_share_bone^2;
% 
% C12(51:200,20:80)=rhob*velocity_longitudinal_bone^2-2*rhob*velocity_share_bone^2;
% %C12=C11-2*C66;
% C21(:,:)=C12(:,:);
% 
% %%
% %中間ゾーン
% lowx(51,20:80)=(rhow+rhob)/2;
% lowx(201,20:80)=(rhow+rhob)/2;
% lowy(51:200,21)=(rhow+rhob)/2;
% lowy(51:200,81)=(rhow+rhob)/2;

inputpointx=1;%音波入力位置　
inputpointy=single(ny/2:ny/2);%音波入力位置

figure;
%% FDTD
 tic;
for s=1:nt
    if s<=n
        Txx(inputpointx,inputpointy)=wave(s);
        Tyy(inputpointx,inputpointy)=wave(s);
    else
        Txx(inputpointx,inputpointy)=0;
        Tyy(inputpointx,inputpointy)=0;
    end
    
%FDTDの式 
        
        Ux(2:nx,1:ny)=Ux(2:nx,1:ny)+dt./lowx(2:nx,1:ny).*((Txx(2:nx,1:ny)-Txx(1:nx-1,1:ny))/dx+(Txy(2:nx,2:ny+1)-Txy(2:nx,1:ny))/dy);
        Uy(1:nx,2:ny)=Uy(1:nx,2:ny)+dt./lowy(1:nx,2:ny).*((Tyy(1:nx,2:ny)-Tyy(1:nx,1:ny-1))/dy+(Txy(2:nx+1,2:ny)-Txy(1:nx,2:ny))/dx);
    
        Txx(1:nx,1:ny)=Txx(1:nx,1:ny)+dt.*(C11(1:nx,1:ny).*(Ux(2:nx+1,1:ny)-Ux(1:nx,1:ny))/dx+C12(1:nx,1:ny).*(Uy(1:nx,2:ny+1)-Uy(1:nx,1:ny))/dy);
        Tyy(1:nx,1:ny)=Tyy(1:nx,1:ny)+dt.*(C21(1:nx,1:ny).*(Ux(2:nx+1,1:ny)-Ux(1:nx,1:ny))/dx+C22(1:nx,1:ny).*(Uy(1:nx,2:ny+1)-Uy(1:nx,1:ny))/dy);
        Txy(2:nx,2:ny)=Txy(2:nx,2:ny)+dt.*C66(2:nx,2:ny).*((Ux(2:nx,2:ny)-Ux(2:nx,1:ny-1))/dy+(Uy(2:nx,2:ny)-Uy(1:nx-1,2:ny))/dx);
        

       
    U=sqrt(Ux.^2+Uy.^2);
    imagesc(U)
%     T=sqrt(Txx.^2+Tyy.^2);
%     imagesc(T)
    colorbar
    axis equal

    M=getframe
    
    
    clc
    s
end
 toc; disp('End');