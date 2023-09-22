%--------------------------------------------------------------------------
% Matlab program  : BIRECTv
% Written by : Naziheddine Belkacem    (naziheddine@yahoo.fr)
% and : Lakhdar Chiter   (lchiter@univ-setif.dz)
% Created on : 0510/2021
% Purpose    : DIRECT optimization algorithm for box constraints.
%--------------------------------------------------------------------------
% clc;close all;
global g example data_ex data_fex f_eval cont
%%%%%%%%%%%%%%%%%%%%%%%  test  example     %%%%%%%%%%%%%%%%%%%%%%%%%%%
% example=9;  % 11 27 39 51
%%%%%%%%%%%%%%%%%%%%%%%  input data   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmaxit=100000;% the maximal number of iterations
nmaxfun=500000;% the maximal number of function evaluations
epsilon=1e-4;% tolerance
pemn=1e-4;%value of pe (percent error)
nafich=9;%(nafich +1) The iteration which you want to display
kdiv=2;%>=2 number of divisions of the graphic windows
graph=0;% 1 for the graphic representation, and 0 if else
%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%
if example==1
% 1 Ackley function (modified) [-15 35]*[-15 35],  f*=0 x*=[0 0]
f_star=0;
%f_star=0.00000837;
%f_star=0.000026;
A0=[-15;-15];
B0=[32.1;32.1];
% A0=[-32.768;-32.768];
% B0=[32.768;32.768];
name='Ackley';
elseif example==2
%   function Ackley dim 5
f_star=0;
% f_star=0.000026;
% f_star=0.0000254
A0=[-15;-15;-15;-15;-15];
B0=[32.1;32.1;32.1;32.1;32.1];
name='Ackley dim 5';
elseif example==3
% 3  function Ackley dim 10
f_star=0;
%f_star=0.000026;
% A0=[-15;-15;-15;-15;-15;-15;-15;-15;-15;-15];
% B0=[32.;32.0;32.2;32.2;32.2;32.2;32.2;32.2;32.2;32.2];
A0=-15*ones(10,1);
B0=32.1*ones(10,1);
name='Ackley dim 10';
elseif example==4
%  Beale function [-4.5 4.5]*[-4.5 4.5], f*=0,x*=[3, 0.5]
f_star=0;
% f_star=0.000096;
A0=[-4.5;-4.5];
B0=[4.5;4.5];
name='Beale';
elseif example==5
%  Bohachevsky1 function (modified) [-100 110] f*=0, x*=[0 0] 
f_star=0;
% f_star=0.0000402;
A0=[-100;-100];
B0=[110.7;110.7];
name='Bohachevsky1';
elseif example==6
%  Bohachevsky2 function (modified) [-100 110] f*=0, x*=[0 0]
f_star=0;
%f_star=0.0000335;
A0=[-100;-100];
B0=[110.7;110.7];
name='Bohachevsky2';
elseif example==7
%  Bohachevsky3 function (modified) [-100 110] f*=0, x*=[0 0]
f_star=0;
%f_star=0.0000367;
A0=[-100;-100];
B0=[110.7;110.7];
name='Bohachevsky3';
elseif example==8
  %  Booth function [-10,10]*[-10,10] f*=0, x*=[1 , 3]
  f_star=0;
% f_star=0.0000610;
  A0=[-10;-10];
  B0=[10.1;10.1];
  name='Booth';
elseif example==9
  %  Branin [-5 10]*[0 15] f*=0.398 ,x*=3vecteurs
  f_star=0.39789;% la valeur minimal th?orique de f
  A0=[-5;0];
  B0=[10;15];
  name='Branin';
elseif example==10
  % Colville function dim 4
   f_star=0.0;
% f_star=0.0000982;
  A0=[-10;-10;-10;-10];
  B0=[10;10;10;10];
  name='Colville';
elseif example==11 
  %  Dixon & Price [-10 10]*[-10 10] f*=0
%   f_star=0.0000484;% la valeur minimal th?orique de f
    f_star=0;
% f_star=0.0000484;
  A0=[-10;-10];
  %B0=[12.0;12.0];
  %B0=[11.1;11.1];
  B0=[10.4554;10.4554];
  name='Dixon & Price dim2';
elseif example==12
    %  Dixon & Price [-10 10]*[-10 10] f*=0
%   f_star=0.0000715;% la valeur minimal th?orique de f
  f_star=0;  
% f_star=0.000072;
%  A0=-10.0*ones(5,1);
%  B0=10.O*ones(5,1);
  A0=-10.40*ones(5,1);
  B0=12.301*ones(5,1);
%  A0=-10.354*ones(5,1);
%   B0=13.10*ones(5,1);
  name='Dixon & Price dim5';
elseif example==13
    %  Dixon & Price [-10 10]*[-10 10] f*=0
%        f_star=0.0000952;% la valeur minimal th?orique de f
  f_star=0;
  A0=-10*ones(10,1);
  B0=12.0*ones(10,1);
  name='Dixon & Price dim10';
elseif example==14
  %  Easom function [-100 100]*[-100 100] f*=-1  x*=[pi pi]
%   f_star=-1;% la valeur minimal th?orique de f
% f_star=-0.99999;
f_star=-1;
  A0=[-100;-100];
  B0=[100;100];
  name='Easom';
elseif example==15
  %  Goldstein and Price [-2 2]*[-2 2] xmin=[0 -1],f*=3 
  f_star=3.0;
%  f_star=3.00019;
  A0=[-2;-2];
  B0=[2;2];
  name='Goldstein';  
elseif example==16
  %  Griewank function [-100 100]*[-100 100] f*=0 x*=[0 0] 
%   f_star=0;
   f_star=0.000000776;
  A0=[-600;-600];
%   B0=[692.2;692.2];
 B0=[700;700];
  name='Griewank';
elseif example==17
 %  Hartman  
  f_star=-3.86278;
  A0=[0;0;0];
  B0=[1;1;1];
  name='Hartman dim 3';
elseif example==18
    %  Hartman  
   f_star=-3.32237;
%  f_star=-3.86242;
% f_star=-3.32206;
  A0=[0;0;0;0;0;0];
  B0=[1;1;1;1;1;1];
  name='Hartman dim 6';
elseif example==19
  %  hump [-5 5]*[-5 5] xmin=[+-0.08984201 -+0.71265640],f*=-1.0316284535
  f_star=-1.03163;
% f_star=-1.03154;
  A0=[-5;-5];
  B0=[5;5];
  name='Hump';  
elseif example==20
  %   Levy function [-10 10]*[-10 10] f*=0 x*=[1 1]
   f_star=0;
% f_star=0.0000909;
  A0=-10.00*ones(2,1);
  B0=10.51*ones(2,1);
  name='Levy 2';
elseif example==21
  %   Levy function [-10 10]*[-10 10] f*=0 x*=[1 1]
  f_star=0;
% f_star=0.0000183;
  A0=-10*ones(5,1);
  B0=10.5*ones(5,1);
  name='Levy 5';
elseif example==22
  %   Levy function [-10 10]*[-10 10] f*=0 x*=[1 1]
   f_star=0;
% f_star=0.0000355;
  A0=-10*ones(10,1);
  B0=10.5*ones(10,1);
  name='Levy 10';    
elseif example==23
  % 14  Matyas function [-10 15]*[-10 15] f*=0 x*=[0 0]
   f_star=0;
% f_star=0.0000271;
  A0=[-10;-10];
  B0=[15;15];
  name='Matyas';
elseif example==24
  % 15  Michalewics function [0 pi]*[0 pi] f*=-1.8013 x*=[2.20 1.57]
  %f_star=-1.80130;
 f_star=-1.80130341009855;
  %f_star=-1.80118;
  A0=[0;0];
  B0=[pi;pi];
  name='Michalewics 2';
elseif example==25
  %   Michalewics 5 function [0 pi]*[0 pi] f*=-1.8013 x*=[2.20 1.57]
    f_star=-4.687658;
%     f_star=-4.687658179088148
    %f_star=-4.68736;
 %  f_star=-4.68721
%   A0=zeros(5,1);
   A0=1.04*ones(5,1);
  B0=pi*ones(5,1);
  name='Michalewics 5';
elseif example==26
  %   Michalewics 5 function [0 pi]*[0 pi] f*=-1.8013 x*=[2.20 1.57]
     f_star=-9.66015;
%       f_star=-8.587;
%     f_star=-7.32661;
%        f_star=-7.87910
%    f_star=-7.84588
  A0=zeros(10,1);
  B0=pi*ones(10,1);
  name='Michalewics 10';
elseif example==27
  %   Perm
   f_star=0;
% f_star=0.00203;
%   A0=-4.5*ones(4,1);
%   B0=4.42*ones(4,1);
A0=-4.0*ones(4,1);
B0=4.0*ones(4,1);
  name='Pern 4';
elseif example==28
  %   Powell 4
  f_star=0;
% f_star=0.0000486;
  A0=-4*ones(4,1);
  B0=5.00*ones(4,1);
  name='Powell 4';
elseif example==29
  %   Powell 8
%  f_star=0.0000971;
 f_star=0;
  A0=-4.00*ones(8,1);
  B0=4.01*ones(8,1);
  name='Powell 8';
elseif example==30
 %   Power sum 4
  f_star=0;
% f_star=0.00009;
    A0=zeros(4,1);
%   A0=1.0*ones(4,1);
%     B0=4.001*ones(4,1);
   B0=4*ones(4,1);
% B0=4.0*ones(4,1);
  name='Power 4';  
elseif example==31
  %   Rastrigin 2 function [-5.12 6.12]*[-5.12 6.12] f*=0 x*=[0 0]
   f_star=0;
%  f_star=0.0000481;
  A0=[-5.12;-5.12];
  B0=[6.12;6.12];
  name='Rastrigin 2'; 
elseif example==32
  %   Rastrigin 5 function [-5.12 6.12]*[-5.12 6.12] f*=0 x*=[0 0]
    f_star=0;
%  f_star=0.0000118;
  A0=-5.12*ones(5,1);
   B0=5.30*ones(5,1);
%  B0=5.12*ones(5,1);
  name='Rastrigin 5'; 
elseif example==33
  %   Rastrigin 10 function [-5.12 6.12]*[-5.12 6.12] f*=0 x*=[0 0]
   f_star=0;
% f_star=0.0000236;
  A0=-5.120*ones(10,1);
  B0=5.1202*ones(10,1);
% A0=-5.12*ones(10,1);
%   B0=5.12*ones(10,1);
  name='Rastrigin 10';   
elseif example==34
  % Rosenbrock 2 function [-5 10]*[-5 10] f*=0  x*=[1 1]
 f_star=0;
% f_star=0.0000965;
  A0=-5*ones(2,1);
  B0=10*ones(2,1);
  name='Rosenbrock 2';
elseif example==35
  % Rosenbrock 5 function [-5 10]*[-5 10] f*=0  x*=[1 1]
   f_star=0;
% f_star=0.0000241;
  A0=-5*ones(5,1);
  B0=10*ones(5,1);
  name='Rosenbrock 5';
elseif example==36
  % Rosenbrock 10 function [-5 10]*[-5 10] f*=0  x*=[1 1]
  f_star=0;
% f_star=0.0000542;
%   A0=-5*ones(10,1);
%   B0=10*ones(10,1);
  A0=-5*ones(10,1);
  B0=10.1*ones(10,1);
  name='Rosenbrock 10'; 
elseif example==37
  %  Schwefel 2 function [-500 500]*[-500 500] f*=0  x*=[0 0]
   f_star=0;
% f_star=0.0000564;
  A0=-519*ones(2,1);
  B0=519*ones(2,1);
  name='Schwefel 2';
elseif example==38
  %  Schwefel 5 function [-500 500]*[-500 500] f*=0  x*=[0 0]
  f_star=0;
% f_star=0.0000641;
  A0=-519*ones(5,1);
  B0=519*ones(5,1);
  name='Schwefel 5';
elseif example==39
  %  Schwefel 10 function [-500 500]*[-500 500] f*=0  x*=[0 0]
%   f_star= 0.00016146;
%      f_star= 0.0000013;
  f_star=0;
  A0=-500*ones(10,1);
  B0=600*ones(10,1);
  name='Schwefel 10';
elseif example==40
  %  Shekel m=5 dim=4
   f_star=-10.15320;
% f_star=-10.15307;
  A0=zeros(4,1);
  B0=10*ones(4,1);
  name='Shenkel m=5 dim=4';
elseif example==41
  %  Shenkel m=7 dim=4
  f_star=-10.40294;
% f_star=-10.40269;
  A0=zeros(4,1);
  B0=10*ones(4,1);
  name='Shenkel m=7 dim=4';
elseif example==42
  %  Shenkel m=10 dim=4
  f_star=-10.53641;
% f_star=-10.53618;
  A0=zeros(4,1);
  B0=10*ones(4,1);
  name='Shenkel m=10 dim=4';
  
elseif example==43
  %shubert [-10 10]*[-10 10] f*=-186.7309
 f_star=-186.73091;
% f_star=-186.72441;
  A0=[-5.12;-5.12];
  B0=[5.12;5.12];
%   A0=[-10;-10];
%   B0=[10;10];
  name='shubert';
elseif example==44
  % Sphere function [-5.12 6.12]*[-5.12 6.12] f*=0  x*=[0 0]
  f_star=0;
% f_star=0.0000115;
  A0=[-5.12;-5.12];
  B0=[6.12;6.12];
%   A0=[-5;-5];
%   B0=[5;5];

  name='Sphere 2';
elseif example==45
  % Sphere function [-5.12 6.12]*[-5.12 6.12] f*=0  x*=[0 0]
   f_star=0;
% f_star=0.0000287;
  A0=-5.12*ones(5,1);
  B0=6.12*ones(5,1);
  name='Sphere 5';
elseif example==46
  % Sphere function [-5.12 6.12]*[-5.12 6.12] f*=0  x*=[0 0]
   f_star=0;
% f_star=0.0000574;
  A0=-5.12*ones(10,1);
  B0=6.12*ones(10,1);
  name='Sphere 10'; 
elseif example==47
    % Sum squares 2
    f_star=0;
% f_star=0.0000794;
  A0=-10*ones(2,1);
  B0=11.5*ones(2,1);
%   A0=-10*ones(2,1);
%   B0=15*ones(2,1);
  name='Sum squares 2';
elseif example==48
    % Sum squares 5
   f_star=0;
%  f_star=0.0000397;
  A0=-10*ones(5,1);
  B0=10.5*ones(5,1);
  name='Sum squares 5';
elseif example==49
    % Sum squares 10
   f_star=0;
% f_star=0.00000911;
  A0=-10*ones(10,1);
  B0=10.5*ones(10,1);
  name='Sum squares 10';
elseif example==50
    % Trid 6
   f_star=-50.0;
% f_star=-49.99512;
  A0=-36.5*ones(6,1);
  B0=36.5*ones(6,1);
  name='Trid 6';
elseif example==51
    % Trid 10
    f_star=-210.0;
% f_star=-209.98007; 
%    A0=-100*ones(10,1);
%    B0=100*ones(10,1);
   A0=-120*ones(10,1);
   B0=120*ones(10,1);
  name='Trid 10'; 
elseif example==52
  % Zakharov function [-5 11]*[-5 11] f*=0  x*=[0 0]
   f_star=0;
% f_star=0.0000288; 
  A0=-5*ones(2,1);
  B0=12*ones(2,1);
  name='Zakharov 2';
elseif example==53
  % Zakharov function [-5 11]*[-5 11] f*=0  x*=[0 0]
  f_star=0;
%  f_star=0.0000644; 
  A0=-5*ones(5,1);
   B0=10.01774*ones(5,1);
% B0=11*ones(5,1);
  name='Zakharov 5';
elseif example==54
  % Zakharov function [-5 11]*[-5 11] f*=0  x*=[0 0]
  f_star=0;
% f_star=9.41133;
  A0=-5*ones(10,1);
  B0=13*ones(10,1);
  name='Zakharov 10';

end
dim=length(A0);
disp(['the name of the chosen test function: ' name]) 
disp(['The lower bound of the domain is: ' num2str(A0') ' The upper bound is: ' num2str(B0')])
disp('********************************************************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cv=diag((B0-A0));
g=@(x)(A0+Cv*x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=zeros(dim,1);B=ones(dim,1);
a=B/3;b=2/3*B;
a=B/3;b=2/3*B;
a=2/3*B;b=B/3;
a=A;b=B;
aa=(b-a);
a=a+aa/3;
b=b;

if dim==2 && graph==1
subplot(kdiv,kdiv,1)
plotrect(A,B,a,b);
end
AA1=[A]; BB1=[B]; aa1=[a]; bb1=[b];
L=[f(g(a))];U=[f(g(b))];nor=norm(g(B)-g(a),2);
[fmin i]=min([L U]);
% a1b1=[a b];nor1=norm(B-a,2);
 a1b1=[a b];nor1=norm(B-a,inf);
xmin=(a1b1(:,i));
D=min([L U]);
Fmin=fmin;
ll=1;
splot=1;
it=1;
XY=[g(a) g(b)];%%%%%%points for figure
f_eval=2;pe=1;
data_ex=[b];
data_fex=[f(g(b))];
while it<=nmaxit && f_eval<=nmaxfun && pe>pemn
%%%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%  %%%%%%%%  %%%%%%  %%%%  %%%%% %%
cont=0;
    %%%
    A=AA1(:,ll);
    B=BB1(:,ll);
    af=aa1(:,ll);
    bf=bb1(:,ll);
    Lf=L(ll);
    Uf=U(ll);
    %%%
    AA1(:,ll)=[];
    BB1(:,ll)=[];
    aa1(:,ll)=[];
    bb1(:,ll)=[];
    L(ll)=[];
    U(ll)=[];nor(ll)=[];D(ll)=[];nor1(ll)=[];
    %%%
    
  for l=1:length(ll) 
  [A1 B1 a1 b1 l1 u1]=Division(A(:,l),B(:,l),af(:,l),bf(:,l),Lf(:,l),Uf(:,l));
  AA1 = [AA1 A1];
  BB1 = [BB1 B1]; 
  aa1 = [aa1 a1]; 
  bb1 =[bb1 b1];

   
L =[L l1];
U =[U u1];
xx=g(B1(:,1))-g(A1(:,1));
nr=norm([xx(1), 2*xx(2)/3]);
[mn ind]= max(B1(:,1)-A1(:,1));

aa=A1(:,1);aa(ind)=aa(ind)+mn/2;
bb=B1(:,1)-(B1(:,1)-A1(:,1))/3;

% nr1=norm(bb-aa);
nr1=2*norm((A1(:,1))-(B1(:,1)))/3;
%nr1=norm((A1(:,1))-(B1(:,1)));
% nr1=2*norm((A1(:,1))-(B1(:,1)),inf)/3;
I=find(abs(nor1-nr1)<=1e-7,1,'first');
 if length(I)==1
     nr1=nor1(I);
 end
XY=[XY,g(a1(:,1)),g(a1(:,2)) g(b1(:,1)),g(b1(:,2))];%%%%%%% points for figure density
nor=[nor nr nr];nor1=[nor1 nr1 nr1];

%%%%%%%%%% min %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a1b1=[a1 b1];
l1u1=[l1 u1];
[l1u1m l1u1i]=min(l1u1);
if l1u1m<=fmin
    xmin=(a1b1(:,l1u1i));
    fmin=l1u1m;
end
Fmin=[Fmin fmin];
D=[D min([l1;u1])];
%%%%%%%%%%%   min  or  max   %%%%%%%%%%%%%%%%%%%%%%%%%%
%f_eval=f_eval+2;
%cont=cont+2;
%%%%%%%%%%% condition  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if f_star==0
       pe=(f(g(xmin)));
    else
       pe=(f(g(xmin))-f_star)/abs(f_star);
    end
% %%%%%%%%%      condition pe<e-4      %%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%
%      if pe<pemn
%     disp(['The value of the percent error is reached: pe= ' num2str(pe) ' at iteration ' num2str(it) ])
%          break
%      end
% %%%%%%%%  condition maximum number of function evaluations %%%%%%%%%%%%%%%% %%%
     if f_eval>nmaxfun
         disp(['The value of maximum number of evaluations is reached: f_eval= ' num2str(f_eval) ' at Iteration ' num2str(it-1) ])
         break
     end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%

  end
  %%%%%%%%  %%%%%%%%%%  %%%%%%%  figures   %%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%
  if dim==2 && graph==1  
  if splot<=kdiv^2-1 
       figure(1)
       subplot(kdiv,kdiv,splot+1)
       for j=1:length(AA1)
        plotrect(AA1(:,j),BB1(:,j),aa1(:,j),bb1(:,j));
       end
       splot=splot+1;
  end
  end
 fprintf('Iteration: %4i   fmin: %15.10f    f evals: %8i\n',it,fmin,cont);
 %%%%%%    Difinir les rectangles optimaux   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  ll=po_opt(nor1,D,fmin,epsilon);%ll=sort(ll,'descend');

%%%%%%%%%%%%%%%%% figure 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dim==2 && graph==1
if it==nafich
figure
  hold on
  for i=1:length(ll)
  A=AA1(:,ll(i));B=BB1(:,ll(i));
  x=[A(1) B(1) B(1) A(1) A(1)];
  y=[A(2) A(2) B(2) B(2) A(2)];
  fill(x,y,'y')
  end
    for j=1:size(AA1,2)%%
    plotrect2(AA1(:,j),BB1(:,j),aa1(:,j),bb1(:,j));
    end 
   title(['Yellow rectangles will be divided next iteration: ' num2str(it+1)])
    figure
    plot(nor1,D,'.r','MarkerSize',10)
    hold on
    plot(nor1(ll),D(ll),'.-b','MarkerSize',12)
    %axis square
    ylabel('min(f(l),f(u))')
    xlabel('Rectangle size')
    legend('Non Potentially Optimal','Potentially Optimal')
    title(['Potentially Optimal & Non Potentially Optimal at Iteration: ' num2str(it+1)])
end
end
it=it+1;
end
disp('******************************* End of Iterations *************************************************')
disp(['after ' num2str(it-1) ' Iterations xmin=[' num2str(g(xmin)','%15.10f') '] , fmin=' num2str(fmin,'%15.10f')])
disp(['f evals = ' num2str(f_eval),' percent error = ' num2str(pe)])
disp('********************************************************************************')
disp(['We recall that : nmaxit= ' num2str(nmaxit),'  nmaxfun= ' num2str(nmaxfun)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dim==2 && graph==1
figure
plot(Fmin,'r')
xlabel('Iterations')
ylabel('f_{min}')
title('Convergence plot')
%%%%%%  %%%%%%%%%  %%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%
figure
h1=max((B0(1)-A0(1))/100,.1);h2=max((B0(2)-A0(2))/100,.1);
[X,Y]= meshgrid(A0(1):h1:B0(1),A0(2):h2:B0(2));
syms x y;
fr=f([x;y]);
disp('the formula of f is: ')
disp(vpa(fr,4))
ff=inline(fr);
Z=ff(X,Y);
meshc(X,Y,Z);
title('Graph of the test function')
figure
hold on
x=[A0(1) B0(1) B0(1) A0(1) A0(1)];
y=[A0(2) A0(2) B0(2) B0(2) A0(2)];
plot(x,y,'-g')
plot(XY(1,:),XY(2,:),'+r','MarkerSize',2)
xx=g(xmin);
plot(xx(1),xx(2),'*b','MarkerSize',3)
title('Scatter plot : points t & v (red color), xmin (blue color)')
end
disp('******************************* end of program *************************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
