%%%
%%% EM Algorithm with Von-mises Distribution by Spturtle (Taejin Park)
%%% This example code demonstrate EM algorithm with Von-mises Distribution.
%%% In this example, 6 sound sources are mixed with arbitary panning
%%% ratio and each TF bins are described with Inter Channel Level Sum Vector (ILVS).
%%% Then, ILVSs are clustered with EM algorithm with Von-mises
%%% Distribution. 


clc;
clear all;
close all;

%%% Bessel function Matrix %%%
xtt=-100:0.01:100;
A=besseli(1, xtt )./besseli( 0, xtt );

%%% Mode Selection %%%
%%% Mode = 1;    5.1ch
%%% Mode = 2;    7.1ch
Mode = 2;

if Mode ==1 
%%% Axis number 
AL=4;    
theta_axis(1)=2*pi*45/360;
theta_axis(2)=2*pi*135/360;
theta_axis(3)=2*pi*225/360;
theta_axis(4)=2*pi*315/360;
end
if Mode ==2 
%%% Axis number 
AL=6;
theta_axis(1)=2*pi*00/360;
theta_axis(2)=2*pi*60/360;
theta_axis(3)=2*pi*120/360;
theta_axis(4)=2*pi*180/360;
theta_axis(5)=2*pi*240/360;
theta_axis(6)=2*pi*300/360;
end

%%% Define Transform Matrix for ILVS
T=zeros(2,2,AL);
d=zeros(5,3);

ang_num=4;
for f=1:AL
T(:,:,f)=[cos(theta_axis(f))  -sin(theta_axis(f)) ; sin(theta_axis(f)) cos(theta_axis(f)) ];    
end

%%% Parameter 
win_size=2048*1;
ap=0.3;
sigma= 0.09;
s = 3; % number of sound
options.subsampling = 1; % Sub Sampling
xf=1:win_size;
wf=sin( (xf-0.5)*pi/ (win_size));

%%% Length of reading buffer
n=48000*8

%%% Read Wavefiles 
[dadda,fs] = wavread('sep_test/bg2.wav',n); 
[jm,fs] = wavread('sep_test/jazz_music.wav',n); 
[bg,fs] = wavread('sep_test/bg2.wav',n); 
[xsrc1(:,1),fs] = wavread('sep_test/BB_bell.wav',n);
[xsrc1(:,2),fs] = wavread('sep_test/es01_mono.wav',n);
[xsrc1(:,3),fs] = wavread('sep_test/car_honk.wav',n);
[xsrc1(:,4),fs] = wavread('sep_test/BB_gus.wav',n);
[xsrc1(:,5),fs] = wavread('sep_test/BB_gun.wav',n);
[xsrc1(:,6),fs] = wavread('sep_test/drum.wav',n);


%%% MIX PRESET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%PRESET 060
if Mode==2;
x1(:,1)=0.5*xsrc1(:,3)+0.7*bg+0.2*xsrc1(:,2);
x1(:,2)=0.5*xsrc1(:,3)+0.7*bg+0.3*xsrc1(:,5);
x1(:,3)=0.2*xsrc1(:,1)+0.7*bg+0.7*xsrc1(:,5);
x1(:,4)=0.8*xsrc1(:,1)+0.7*bg+0.2*xsrc1(:,6)+0.3*xsrc1(:,5);
x1(:,5)=0.5*xsrc1(:,4)+0.7*bg+0.8*xsrc1(:,6);
x1(:,6)=0.5*xsrc1(:,4)+0.7*bg+0.8*xsrc1(:,2);
end

xz=(  (xsrc1(:,1)+xsrc1(:,3)+xsrc1(:,4)+xsrc1(:,5))/4 );

%%% MIX PRESET END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FN=2*floor( n/win_size )-1;
color=4;

gain_mat=zeros(5,3);
delay_mat=zeros(5,3);
overall_gain=3;

mdct_mic=zeros(win_size/2,FN,AL);
mdct_mic_p=zeros(win_size/2,FN,AL);
mdct_xz=zeros(win_size/2,FN,1);
mdct_mic_y=zeros(win_size/2,FN);
mdct_mic_x=zeros(win_size/2,FN);
mdct_mic_sum=zeros(win_size/2,FN,1);
r_value=zeros(win_size/2,FN,1);

figure;
theta  = linspace(0,2*pi,100);
r      = sin(2*theta) .* cos(2*theta);
%%% Set Range of Polar Plot

d_gain=200;
r_max  = 4*d_gain;
h_fake = polar(theta,r_max*ones(size(theta)));
hold on;
set(h_fake, 'Visible', 'Off');

polar([theta_axis(1),theta_axis(1)],[0,r_max])
polar([theta_axis(2),theta_axis(2)],[0,r_max])
polar([theta_axis(3),theta_axis(3)],[0,r_max])
polar([theta_axis(4),theta_axis(4)],[0,r_max])
polar([theta_axis(5),theta_axis(5)],[0,r_max])
polar([theta_axis(6),theta_axis(6)],[0,r_max])

THETA=zeros(win_size/2,FN);
RHO=zeros(win_size/2,1);

ct=1;
for r=1:FN;
  r;
for u=1:AL
    u;
mdct_mic(:,r,u)=mdct4( wf(:).*x1(   ( win_size/2*(r-1) +1) : ( win_size/2*(r-1) +win_size) ,u) );
mdct_mic_p(:,r,u)=mdct4( wf(:).*x1(   ( win_size/2*(r-1) +1) : ( win_size/2*(r-1) +win_size) ,u) );
mdct_mic_sum(:,r)= mdct_mic_sum(:,r) +abs( mdct_mic(:,r,u) );
end

mdct_xz(:,r,1)=mdct4( wf(:).*xz(   ( win_size/2*(r-1) +1) : ( win_size/2*(r-1) +win_size) ,1) );
for u=1:AL
mdct_mic_x(:,r)=mdct_mic_x(:,r)+T(1,1,u)*abs(mdct_mic(:,r,u));
mdct_mic_y(:,r)=mdct_mic_y(:,r)+T(2,1,u)*abs(mdct_mic(:,r,u));
end

% plot(mdct_mic_x(:,r),mdct_mic_y(:,r),'.')
[THETA(:,r),RHO(:,r)] = cart2pol(mdct_mic_x(:,r),mdct_mic_y(:,r));
polar(THETA(:,r),d_gain*RHO(:,r),'.');

r_value(:,r)=( RHO(:,r)  )./ ( mdct_mic_sum(:,r) );

for q=1:(win_size/2)
if RHO(q,r)  > 0.01;

r_den(ct)=( RHO(q,r)  )/(mdct_mic_sum(q,r) );

ct=ct+1;

end
end

title('Inter-channel Level Vector Sum','FontSize',15)

end

TH_all=reshape(THETA,FN*win_size/2,1);
RH_all=reshape(RHO,FN*win_size/2,1);
RV_all=reshape(r_value,FN*win_size/2,1);
RV_all=r_den';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cl=8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% X=ang_cum;
alpha=1*ones(cl,1);
c=0.001*ones(cl,1);
%%%% eps is VERY IMPORTANT value ( Have influence on the EM fitting
%%%% performance)

rx=0:0.01:(max(RV_all));

for p=1:AL;
x_res(:,p)=x1(:,p);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VON MISES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
get_cl=6;
mdct_mic_y=zeros(win_size/2,FN);
mdct_mic_x=zeros(win_size/2,FN);
mdct_mic=zeros(win_size/2,FN,AL);
mdct_xz=zeros(win_size/2,FN,1);

ct=1;
for r=1:FN;
  r
for u=1:AL
mdct_mic(:,r,u)=mdct4( wf(:).*x_res(   ( win_size/2*(r-1) +1) : ( win_size/2*(r-1) +win_size) ,u) );
mdct_mic_sum(:,r)= mdct_mic_sum(:,r) +abs( mdct_mic(:,r,u) );
end

mdct_xz(:,r,1)=mdct4( wf(:).*xz(   ( win_size/2*(r-1) +1) : ( win_size/2*(r-1) +win_size) ,1) );
for u=1:AL
mdct_mic_x(:,r)=mdct_mic_x(:,r)+T(1,1,u)*abs(mdct_mic(:,r,u));
mdct_mic_y(:,r)=mdct_mic_y(:,r)+T(2,1,u)*abs(mdct_mic(:,r,u));
end

% plot(mdct_mic_x(:,r),mdct_mic_y(:,r),'.')
[THETA(:,r),RHO(:,r)] = cart2pol(mdct_mic_x(:,r),mdct_mic_y(:,r));
% polar(THETA(:,r),RHO(:,r),'.');

r_value(:,r)=( RHO(:,r)  )./ ( mdct_mic_sum(:,r) );

for q=1:(win_size/2)
if RHO(q,r)  > 0.05;

th_den(ct)=THETA(q,r);

ct=ct+1;

end
end
%title('MDCT 2ch plot')
hold on

end

TH_all=reshape(THETA,FN*win_size/2,1);
TH_all=th_den';
RH_all=reshape(RHO,FN*win_size/2,1);
RV_all=reshape(r_value,FN*win_size/2,1);
RV_all=r_den';

MS=zeros(1,cl);

alpha=1*ones(cl,1);
c=20*ones(1,cl);
%%%% eps is VERY IMPORTANT value ( Have influence on the EM fitting
%%%% performance)
eps=pi/2;

x= -(pi):0.01:(pi);
% figure;
[x3out,y]=hist(TH_all,x);
y3=x3out;

x= -(pi):0.01:(pi);
theta  = linspace(0,2*pi,100);
r      = sin(2*theta) .* cos(2*theta);
%%% Set Range of Polar Plot
r_max  = max((y3));
h_fake = polar(theta,r_max*ones(size(theta)));
hold on;
set(h_fake, 'Visible', 'Off');
hold on
polar(x,log2(y3));
hold on
%%% Initialize MEAN value of Laplacian cluster

MV=linspace(-pi,pi,cl);
m=10*ones(1,cl);

scale=r_max*20;
k=cl;
X=TH_all;

xt=x;
old_MV=0*MV;
ct=0;


while sum(abs(MV-old_MV)) >0.00005  && ct <80
clc;
m;
old_MV=MV
ct=ct+1

%%%%%%%% EM algorithm with Von-mises Distribution 
[E]=exp_step_P(xt,X,cl,alpha,MV,m);
[alpha,MV,m]=max_step_P(xt,X,cl,E,MV,eps,A);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:cl

dXM = x-MV(j);
y= (   1/( sum ( exp(m(j)*cos(xt) ) ) )*exp(m(j)*cos(dXM)  )    );
if j==1
polar(x,scale*y,'m-')
hold on
end

if j==2
polar(x,scale*y,'k-')
hold on
end

if j==3
polar(x,scale*y,'g-')
hold on
end

if j==4
polar(x,scale*y,'c-')
hold on
end

if j==5
polar(x,scale*y,'r-')
hold on
end

if j==6
polar(x,scale*y,'y-')
hold on
end

if j==7
polar(x,scale*y,'m-.')
hold on
end

if j==8
polar(x,scale*y,'k-.')
hold on
end
if j==9
polar(x,scale*y,'g-.')
hold on
end

if j==10
polar(x,scale*y,'c-.')
hold on
end

end
pause(1e-20)
end


for i=1:cl
if m(i) < 0
    MS(i)=MV(i)+pi; 
    if MS(i) > pi
        MS(i) = MS(i)-2*pi;
    end
end
if m(i) > 0
  MS(i)=MV(i);   
end
end

%%% Final graph display

hold off

figure;
h_fake = polar(theta,r_max*ones(size(theta)));
hold on;
set(h_fake, 'Visible', 'Off');
hold on
polar(x,y3);
hold on

for j=1:cl
dXM = x-MV(j);

br(:,j)= (   1/( sum ( exp(m(j)*cos(xt) ) ) )*exp(m(j)*cos(dXM)  )    );
if j==1
polar(x,scale*br(:,j)','m-')
hold on
end

if j==2
polar(x,scale*br(:,j)','k-')
hold on
end

if j==3
polar(x,scale*br(:,j)','g-')
hold on
end

if j==4
polar(x,scale*br(:,j)','c-')
hold on
end

if j==5
polar(x,scale*br(:,j)','r-')
hold on
end

if j==6
polar(x,scale*br(:,j)','y-')
hold on
end

if j==7
polar(x,scale*br(:,j)','m-.')
hold on
end

if j==8
polar(x,scale*br(:,j)','k-.')
hold on
end

if j==9
polar(x,scale*br(:,j)','g-.')
hold on
end

if j==10
polar(x,scale*br(:,j)','c-.')
hold on
end

end

for yt=1:cl;
tip_cl_von(yt)=max( br(:,yt) );
end

[cit,index_von_get]=sort(tip_cl_von,'descend');

pr=index_von_get(1:get_cl);
for tt=1:get_cl
get_br(:,(tt) )= (   1./( sum ( exp(m(pr(tt))*cos(xt) ) ) ).*exp(m(pr(tt))*cos(x-MV(pr(tt)))  )    );
end


x_len=length(x);
pri=zeros(x_len,1);
for tr=1:x_len
    
[ct,ir]=sort(get_br(tr,:),'descend');
    pri(tr)=pr(ir(1));
end




pri_size=size(pri);
MS
sep_out=zeros(win_size/2,FN,cl);
sep_out_div=zeros(win_size/2,FN,cl);
sep_file=zeros(length(x1(:,1)), cl );
ths=linspace(-pi,pi,pri_size(1));


angles=zeros(ang_num,1); 

if Mode==1;
angles(1)=pi/4;
angles(2)=3*pi/4;
angles(3)=-3*pi/4;
angles(4)=-pi/4;;
end 

if Mode==2;
angles(1)=0;
angles(2)=pi/3;
angles(3)=2*pi/3;
angles(4)=-pi;
angles(5)=-2*pi/3;;
angles(6)=-pi/3;;
end 

mdct_write(:,:,1)=mdct_mic(:,:,1);
mdct_write(:,:,2)=mdct_mic(:,:,2);
mdct_write(:,:,3)=mdct_mic(:,:,3);
mdct_write(:,:,4)=mdct_mic(:,:,4);
if Mode==2;
mdct_write(:,:,5)=mdct_mic(:,:,5);
mdct_write(:,:,6)=mdct_mic(:,:,6);
end
for r=1:FN    
    r

[THETA(:,r),RHO(:,r)] = cart2pol(mdct_mic_x(:,r),mdct_mic_y(:,r));
for f=1:win_size/2


fv=ths-THETA(f,r);
[c1,minfv_id]=min(abs(fv));
slotfv_id=pri(minfv_id);
% end

av=angles-THETA(f,r);
[c2,minav_id]=sort(abs(av));

if c2(1) < 0.25
if  r_value(f,r) > 0.8 && r(f,r) < 1.2
    sep_out(f,r,slotfv_id)=mdct_write(f,r,minav_id(1));
else
    sep_out(f,r,slotfv_id)=mdct_write(f,r,minav_id(1))+mdct_write(f,r,minav_id(2))+mdct_write(f,r,minav_id(3));
end

else
sep_out(f,r,slotfv_id)=mdct_write(f,r,minav_id(1))+mdct_write(f,r,minav_id(2));
end

   

end
for u=1:cl
sep_file( ((win_size/2)*(r-1)+1 ) : ((win_size/2)*(r-1)+win_size ) ,u)=sep_file( ((win_size/2)*(r-1)+1 ) : ((win_size/2)*(r-1)+win_size ),u )   +  wf(:).* imdct4( sep_out(:,r,u) );
end

end

for CN=1:cl
 wavwrite(sep_file(:,CN),fs,strcat('sep_test/auto_sep_out',num2str(CN),'.wav'))
end
hold off

x= -(pi):0.01:(pi);
r_max=1000;

scale=8.1152e+03;
scale=3*scale;
figure;
%%%%
p=polar(x,100*log2(y3));
set(p,'linewidth',1.5)
% plot(x,y3);
hold on
p1=polar([theta_axis(1),theta_axis(1)],[0,r_max])
set(p1,'linewidth',1)
p2=polar([theta_axis(2),theta_axis(2)],[0,r_max])
set(p2,'linewidth',1)
p3=polar([theta_axis(3),theta_axis(3)],[0,r_max])
set(p3,'linewidth',1)
p4=polar([theta_axis(4),theta_axis(4)],[0,r_max])
set(p4,'linewidth',1)
p5=polar([theta_axis(5),theta_axis(5)],[0,r_max])
set(p5,'linewidth',1)
p6=polar([theta_axis(6),theta_axis(6)],[0,r_max])
set(p6,'linewidth',1)

hold on
for rr=1:get_cl
j=pr(rr);
% m=r0/(c(j)^2);
dXM = x-MV(j);
% y= exp(m)/(2*pi*sqrt(2*pi*m))*exp(m*cos(dXM));
br(:,j)= (   1/( sum ( exp(m(j)*cos(xt) ) ) )*exp(m(j)*cos(dXM)  )    );
if j==1
h1=polar(x,scale*br(:,j)','m--')
set(h1,'linewidth',1.5)
hold on
end

if j==2
h2=polar(x,scale*br(:,j)','k--')
set(h2,'linewidth',1.5)
hold on
end

if j==3
h3=polar(x,scale*br(:,j)','g--')
set(h3,'linewidth',1.5)
hold on
end

if j==4
h4=polar(x,scale*br(:,j)','c--')
set(h4,'linewidth',1.5)
hold on
end

if j==5
h5=polar(x,scale*br(:,j)','r--')
set(h5,'linewidth',1.5)
hold on
end

if j==6
h6=polar(x,scale*br(:,j)','r--')
set(h6,'linewidth',1.5)
hold on
end

if j==7
h7=polar(x,scale*br(:,j)','m-.')
set(h7,'linewidth',1.5)
hold on
end

if j==8
h8=polar(x,scale*br(:,j)','k-.')
set(h8,'linewidth',1.5)

hold on
end

if j==9
h9=polar(x,scale*br(:,j)','g-.')
set(h9,'linewidth',1.5)
hold on
end

if j==10
h10=polar(x,scale*br(:,j)','c-.')
set(h10,'linewidth',1.5)
hold on
end

end
