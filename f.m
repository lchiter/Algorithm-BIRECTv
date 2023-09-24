% MATLAB coding by: Naziheddine Belkacem & Lakhdar Chiter
% A collection of 54 global optimization test problems from :
% Source:
%  - http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page295.htm
%--------------------------
function y=f(x)
global example
if example ==1
    % Ackley function.
n = 2;
a = 20; b = 0.2; c = 2*pi;
s1 = 0; s2 = 0;
for i=1:n;
   s1 = s1+x(i)^2;
   s2 = s2+cos(c*x(i));
end
y = -a*exp(-b*(1/n*s1)^(1/2))-exp(1/n*s2)+a+exp(1);
%---------------------------
elseif example ==2
    % Ackley function.
    n = 5;
a = 20; b = 0.2; c = 2*pi;
s1 = 0; s2 = 0;
for i=1:n;
   s1 = s1+x(i)^2;
   s2 = s2+cos(c*x(i));
end
y = -a*exp(-b*sqrt(1/n*s1))-exp(1/n*s2)+a+exp(1);
%---------------------------
elseif example==3
    % Ackley function.
    n = 10;
a = 20; b = 0.2; c = 2*pi;
s1 = 0; s2 = 0;
for i=1:n;
   s1 = s1+x(i)^2;
   s2 = s2+cos(c*x(i));
end
% y = -a*exp(-b*sqrt(1/n*s1))-exp(1/n*s2)+a+exp(1);
y = -a*exp(-b*(1/n*s1)^(1/2))-exp(1/n*s2)+a+exp(1);
%---------------------------
elseif example==4
% Beale function.
% The number of variables n = 2.
% 
y = (1.5-x(1)*(1-x(2)))^2+(2.25-x(1)*(1-x(2)^2))^2+(2.625-x(1)*(1-x(2)^3))^2;
%-----------------------------
elseif example==5
% Bohachecsky function 1
y=x(1)^2+2*x(2)^2-0.3*cos(3*pi*x(1))-0.4*cos(4*pi*x(2))+0.7;
%------------------------------
elseif example==6
% Bohachecsky function 2
% The number of variables n = 2. 
y = x(1)^2+2*x(2)^2-0.3*cos(3*pi*x(1))*cos(4*pi*x(2))+0.3;
elseif example==7
% Bohachecsky function 3
% The number of variables n = 2. 
y = x(1)^2+2*x(2)^2-0.3*cos(3*pi*x(1)+4*pi*x(2))+0.3;
%-----------------------------
elseif example==8
% Booth function 
% The number of variables n = 2.
y  = (x(1)+2*x(2)-7)^2+(2*x(1)+x(2)-5)^2;
%---------------------------------------
elseif example==9
% Branin function 
% The number of variables n = 2.
y = (x(2)-(5.1/(4*pi^2))*x(1)^2+5*x(1)/pi-6)^2+10*(1-1/(8*pi))*cos(x(1))+10;
%----------------------------
elseif example==10
% Colville function 
% The number of variables n = 4. 
y  = 100*(x(1)^2-x(2))^2+(x(1)-1)^2+(x(3)-1)^2+90*(x(3)^2-x(4))^2+10.1*((x(2)-1)^2+(x(4)-1)^2)+19.8*(x(2)-1)*(x(4)-1);
%---------------------------------
elseif example==11
% Dixon and Price function.
% The number of variables n should be adjusted below.
% The default value of n = 2, 5, 10.
s1 = 0;n=2;
for j = 2:n;
    s1 = s1+j*(2*x(j)^2-x(j-1))^2;    
end
y = s1+(x(1)-1)^2;
%----------------------------------
elseif example==12
    % Dixon and Price function.
% The number of variables n should be adjusted below.
% The default value of n = 2, 5, 10.
n=5;
term1=(x(1)-1)^2;
sum = 0;
for j = 2:n;
    xold = x(j-1);
    new = j*(2*x(j)^2-xold)^2;
    sum = sum+new;    
end
y = term1+sum;
%----------------------------------
elseif example==13
    % Dixon and Price function.
% The number of variables n should be adjusted below.
% The default value of n = 2, 5, 10.
n=10;
% term1=(x(1)-1)^2;
s1= 0;
for j = 2:n
%     xold = x(j-1);
%     new = j*(2*x(j)^2-xold)^2;
      s1 = s1 + j*(2*x(j)^2 - x(j - 1))^2; 
end
y = s1 + (x(1) - 1)^2;
%----------------------------------
elseif example==14
% Easom function 
% The number of variables n = 2. 
y = -cos(x(1))*cos(x(2))*exp(-(x(1)-pi)^2-(x(2)-pi)^2);
%-----------------------------------
elseif example==15
% Goldstein and Price function 
% The number of variables n = 2. 
a = 1+(x(1)+x(2)+1)^2*(19-14*x(1)+3*x(1)^2-14*x(2)+6*x(1)*x(2)+3*x(2)^2);
b = 30+(2*x(1)-3*x(2))^2*(18-32*x(1)+12*x(1)^2+48*x(2)-36*x(1)*x(2)+27*x(2)^2);
y = a*b;
%---------------------------
elseif example==16
% Griewank function
% The number of variables n should be adjusted below.
% The default value of n =2. 
n = 2;
fr = 4000;
s = 0;
p = 1;
for j = 1:n 
    s = s+x(j)^2; 
end
for j = 1:n
    p = p*cos(x(j)/sqrt(j)); 
end
y = s/fr - p + 1;
%-----------------------------------------
elseif example==17
% Hartmann function 
% The number of variables n = 3. 
a(:,2)=10.0*ones(4,1);
for j=1:2;
   a(2*j-1,1)=3.0; a(2*j,1)=0.1; 
   a(2*j-1,3)=30.0; a(2*j,3)=35.0; 
end
c(1)=1.0;c(2)=1.2;c(3)=3.0;c(4)=3.2;
p(1,1)=0.36890;p(1,2)=0.11700;p(1,3)=0.26730;
p(2,1)=0.46990;p(2,2)=0.43870;p(2,3)=0.74700;
p(3,1)=0.10910;p(3,2)=0.87320;p(3,3)=0.55470;
p(4,1)=0.03815;p(4,2)=0.57430;p(4,3)=0.88280;
sum = 0;
for i=1:4;
   sm=0;
   for j=1:3;
      sm=sm+a(i,j)*(x(j)-p(i,j))^2;
   end
   sum=sum+c(i)*exp(-sm);
end
y = -sum;
%---------------------------
elseif example==18 
% Hartmann function 
% Matlab Code by A. Hedar (Sep. 29, 2005).
% The number of variables n = 6. 
a(1,1)=10.0;a(1,2)=3.0;	a(1,3)=17.0;a(1,4)=3.5;	a(1,5)=1.7;	a(1,6)=8.0;
a(2,1)=0.05;a(2,2)=10.0;a(2,3)=17.0;a(2,4)=0.1;	a(2,5)=8.0;	a(2,6)=14.0;
a(3,1)=3.0;	a(3,2)=3.5;a(3,3)=1.7;	a(3,4)=10.0;a(3,5)=17.0;a(3,6)=8.0;
a(4,1)=17.0;a(4,2)=8.0;a(4,3)=0.05;	a(4,4)=10.0;a(4,5)=0.1;	a(4,6)=14.0;
c(1)=1.0;c(2)=1.2;c(3)=3.0;c(4)=3.2;
p(1,1)=0.1312;p(1,2)=0.1696;p(1,3)=0.5569;p(1,4)=0.0124;p(1,5)=0.8283;p(1,6)=0.5886;
p(2,1)=0.2329;p(2,2)=0.4135;p(2,3)=0.8307;p(2,4)=0.3736;p(2,5)=0.1004;p(2,6)=0.9991;
p(3,1)=0.2348;p(3,2)=0.1451;p(3,3)=0.3522;p(3,4)=0.2883;p(3,5)=0.3047;p(3,6)=0.6650;
p(4,1)=0.4047;p(4,2)=0.8828;p(4,3)=0.8732;p(4,4)=0.5743;p(4,5)=0.1091;p(4,6)=0.0381;
sum = 0;
for i=1:4;
   sm=0;
   for j=1:6;
      sm=sm+a(i,j)*(x(j)-p(i,j))^2;
   end
   sum=sum+c(i)*exp(-sm);
end
y = -sum;
%---------------------------
elseif example==19
% Hump function 
% The number of variables n = 2.
% y=1.0316285+4*x(1)^2-2.1*x(1)^4+x(1)^6/3+x(1)*x(2)-4*x(2)^2+4*x(2)^4;
y = (4-2.1*x(1)^2+(x(1)^4)/3)*x(1)^2+x(1)*x(2)+(-4+4*x(2)^2)*x(2)^2;
%----------------------------------
elseif example==20
% Levy function 
% The number of variables n should be adjusted below.
% The default value of n =2.
% 
n = 2;% 5, 10;
for i = 1:n; z(i) = 1+(x(i)-1)/4; end
sum = sin(pi*z(1))^2;
for i = 1:n-1
    sum = sum+(z(i)-1)^2*(1+10*(sin(pi*z(i)+1))^2);
end 
y = sum+(z(n)-1)^2*(1+(sin(2*pi*z(n)))^2);
%----------------------------------------------
elseif example==21
    % Levy function 
% The number of variables n should be adjusted below.
% The default value of n =2.
% 
n = 5;% 5, 10;
for i = 1:n; 
 z(i) = 1+(x(i)-1)/4; 
end
s = sin(pi*z(1))^2;
for i = 1:n-1
    s = s+(z(i)-1)^2*(1+10*(sin(pi*z(i)+1))^2);
end 
y = s+(z(n)-1)^2*(1+(sin(2*pi*z(n)))^2);
elseif example==22
    % Levy function 
% The number of variables n should be adjusted below.
% The default value of n =2.
% 
n = 10;% 5, 10;
for i = 1:n; z(i) = 1+(x(i)-1)/4; end
sum = sin(pi*z(1))^2;
for i = 1:n-1
    sum = sum+(z(i)-1)^2*(1+10*(sin(pi*z(i)+1))^2);
end 
y = sum+(z(n)-1)^2*(1+(sin(2*pi*z(n)))^2);
elseif example==23
% Matyas function 
% The number of variables n =2.
% 
y = 0.26*(x(1)^2+x(2)^2)-0.48*x(1)*x(2);
%-------------------------------------------
elseif example==24
% Michalewicz function 
% The number of variables n should be adjusted below.
% The default value of n =2.
% 
n = 2; 
m = 10;
sum = 0;
for i = 1:n;
    sum = sum+sin(x(i))*(sin(i*x(i)^2/pi))^(2*m);
end
y = -sum;
%-----------------------------------------
elseif example==25
 % Michalewicz function 
% The number of variables n should be adjusted below.
% The default value of n =2.
% 
m = 10;
n = 5; 
s = 0;
for i = 1:n
    new = sin(x(i))*(sin(i*x(i)^2/pi))^(2*m);
    s = s + new;
end
y = -s;
%----------------------------------------------- 
elseif example==26
    % Michalewicz function 
% The number of variables n should be adjusted below.
% The default value of n =2.
% 
n = 10; 
m = 10;
s = 0;
for i = 1:n;
    s = s+sin(x(i))*(sin(i*x(i)^2/pi))^(2*m);
end
y = -s;
%-----------------------------------------------    
elseif example==27
% Perm function 
% The number of variables n should be adjusted below.
% The default value of n = 4.
% 
n = 4;
% b=10
b = 0.5;
outer = 0;
for k = 1:n
    inner = 0;
    for j = 1:n
        inner = inner+(j^k+b)*((x(j)/j)^k-1);
%         inner = inner + (j + b)*(x(j)^k - (1/j)^k);
    end
    outer = outer+inner^2;
end
y = outer;
%------------------------------------------------
elseif example==28
% Powell function 
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% 
n = 4;%, 8;
m = n;
% for i = 1:m/4
%     fvec(4*i-3) = x(4*i-3)+10*(x(4*i-2));
%     fvec(4*i-2) = sqrt(5)*(x(4*i-1)-x(4*i));
%     fvec(4*i-1) = (x(4*i-2)-2*(x(4*i-1)))^2;
%     fvec(4*i)   = sqrt(10)*(x(4*i-3)-x(4*i))^2;
sum = 0;
for i = 1:m/4
 %for i = 1:2
    term1 = (x(4*i-3)+10*x(4*i-2))^2;
	term2 = 5*(x(4*i-1)-x(4*i))^2;
	term3 = (x(4*i-2)-2*x(4*i-1))^4;
	term4 = 10*(x(4*i-3)-x(4*i))^4;
	sum = sum + term1 + term2 + term3 + term4;
end;
% fvec = fvec';
% y = norm(fvec)^2;
y = sum;
%----------------------------
elseif example==29
% Powell function 
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% 
  n = 8;%, 8;
  m = n;
sum = 0;
for i = 1:m/4
 %for i = 1:2
    term1 = (x(4*i-3)+10*x(4*i-2))^2;
	term2 = 5*(x(4*i-1)-x(4*i))^2;
	term3 = (x(4*i-2)-2*x(4*i-1))^4;
	term4 = 10*(x(4*i-3)-x(4*i))^4;
	sum = sum + term1 + term2 + term3 + term4;
%      fvec(4*i-3) = x(4*i-3)+10*(x(4*i-2));
%      fvec(4*i-2) = sqrt(5)*(x(4*i-1)-x(4*i));
%      fvec(4*i-1) = (x(4*i-2)-2*(x(4*i-1)))^2;
%      fvec(4*i)   = sqrt(10)*(x(4*i-3)-x(4*i))^2;
end;
%   fvec = fvec';
%   y = norm(fvec)^2;
y = sum;
%----------------------------
elseif example==30
% Power Sum function 
% The number of variables n should be adjusted below.
% The default value of n = 4.
% 
n = 4;
 b = [8,18,44,114];
s_out = 0;
for k = 1:n
    s_in = 0;
    for j = 1:n
        s_in = s_in+x(j)^k;
    end
    s_out = s_out+(s_in-b(k))^2;
end
y = s_out;
%----------------------------
elseif example==31
% Rastrigin function
% The number of variables n should be adjusted below.
% The default value of n = 2. (2, 5, et 10)
% 
n = 2;% 5, 10;
s = 0;
for j = 1:n
    s = s+(x(j)^2-10*cos(2*pi*x(j))); 
end
y = 10*n+s;
%-------------------------
elseif example==32
% Rastrigin function
% The number of variables n should be adjusted below.
% The default value of n = 2. (2, 5, et 10)
% 
n = 5;% 5, 10;
s = 0;
for j = 1:n
    s = s+(x(j)^2-10*cos(2*pi*x(j))); 
end
y = 10*n+s;
%-------------------------
elseif example==33
% Rastrigin function
% The number of variables n should be adjusted below.
% The default value of n = 2. (2, 5, et 10)
% 
n = 10;% 5, 10;
s = 0;
for j = 1:n
    s = s+(x(j)^2-10*cos(2*pi*x(j))); 
end
y = 10*n+s;
%-------------------------
elseif example==34
% Rosenbrock function
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 2.
% 
n = 2;%, 5, 10;
sum = 0;
for j = 1:n-1;
    sum = sum+100*(x(j)^2-x(j+1))^2+(x(j)-1)^2;
end
y = sum;
%-------------------------------------------
elseif example==35
% Rosenbrock function
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 2.
% 
n = 5;%, 5, 10;
sum = 0;
for j = 1:n-1;
    sum = sum+100*(x(j)^2-x(j+1))^2+(x(j)-1)^2;
end
y = sum;
%-------------------------------------------
elseif example==36
% Rosenbrock function
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 2.
% 
n = 10;%, 5, 10;
sum = 0;
for j = 1:n-1;
    sum = sum+100*(x(j)^2-x(j+1))^2+(x(j)-1)^2;
end
y = sum;
%-------------------------------------------
elseif example==37
% Schwefel function
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 2.
% 
n = 2;%, 5, 10;
s = 0;
for j = 1:n;
   s = s + x(j)*sin(sqrt(abs(x(j))));
end
y = 418.9828872724336*n - s;
%------------------------
elseif example==38
% Schwefel function
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 2.
% 
n = 5;%, 5, 10;
s = 0;
for j = 1:n;
   s = s + x(j)*sin(sqrt(abs(x(j))));
end
y = 418.9828872724336*n - s;
%------------------------
elseif example==39
% Schwefel function
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 2.
% 
   n = 10;%, 5, 10;
   s= 0;
  for j = 1:n;
     s = s + x(j)*sin(sqrt(abs(x(j))));
%       s =s -x(j)*sin(sqrt(abs(x(j))));
 end
%     y = 418.9829*n-s;
%   y = 418.9829*n-sum;
 %   s = sum(-x.*sin(sqrt(abs(x))));
% s = sum(-x(j).*sin(sqrt(abs(x(j)))));
 y = 418.9828872724336*n - s;
%------------------------
elseif example==40
% Shekel function
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n = 4
% The parameter m should be adjusted m = 5,7,10.
% The default value of m = 10.
% 
m = 5;
a = ones(10,4);
a(1,:) = 4.0*a(1,:);
a(2,:) = 1.0*a(2,:);
a(3,:) = 8.0*a(3,:);
a(4,:) = 6.0*a(4,:);
for j = 1:2;
   a(5,2*j-1) = 3.0; a(5,2*j) = 7.0; 
   a(6,2*j-1) = 2.0; a(6,2*j) = 9.0; 
   a(7,j)     = 5.0; a(7,j+2) = 3.0;
   a(8,2*j-1) = 8.0; a(8,2*j) = 1.0;
   a(9,2*j-1) = 6.0; a(9,2*j) = 2.0;
   a(10,2*j-1)= 7.0; a(10,2*j)= 3.6;
end
c(1) = 0.1; c(2) = 0.2; c(3) = 0.2; c(4) = 0.4; c(5) = 0.4;
c(6) = 0.6; c(7) = 0.3; c(8) = 0.7; c(9) = 0.5; c(10)= 0.5;
sum = 0;
for j = 1:m;
   p = 0;
   for i = 1:4
      p = p+(x(i)-a(j,i))^2;
   end
   sum = sum+1/(p+c(j));
end
y = -sum;
%--------------------------------
elseif example==41
% Shekel function
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n = 4
% The parameter m should be adjusted m = 5,7,10.
% The default value of m = 10.
% 
m = 7;
a = ones(10,4);
a(1,:) = 4.0*a(1,:);
a(2,:) = 1.0*a(2,:);
a(3,:) = 8.0*a(3,:);
a(4,:) = 6.0*a(4,:);
for j = 1:2;
   a(5,2*j-1) = 3.0; a(5,2*j) = 7.0; 
   a(6,2*j-1) = 2.0; a(6,2*j) = 9.0; 
   a(7,j)     = 5.0; a(7,j+2) = 3.0;
   a(8,2*j-1) = 8.0; a(8,2*j) = 1.0;
   a(9,2*j-1) = 6.0; a(9,2*j) = 2.0;
   a(10,2*j-1)= 7.0; a(10,2*j)= 3.6;
end
c(1) = 0.1; c(2) = 0.2; c(3) = 0.2; c(4) = 0.4; c(5) = 0.4;
c(6) = 0.6; c(7) = 0.3; c(8) = 0.7; c(9) = 0.5; c(10)= 0.5;
sum = 0;
for j = 1:m;
   p = 0;
   for i = 1:4
      p = p+(x(i)-a(j,i))^2;
   end
   sum = sum+1/(p+c(j));
end
y = -sum;
%--------------------------------
elseif example==42
% Shekel function
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n = 4
% The parameter m should be adjusted m = 5,7,10.
% The default value of m = 10.
% 
m = 10;
a = ones(10,4);
a(1,:) = 4.0*a(1,:);
a(2,:) = 1.0*a(2,:);
a(3,:) = 8.0*a(3,:);
a(4,:) = 6.0*a(4,:);
for j = 1:2;
   a(5,2*j-1) = 3.0; a(5,2*j) = 7.0; 
   a(6,2*j-1) = 2.0; a(6,2*j) = 9.0; 
   a(7,j)     = 5.0; a(7,j+2) = 3.0;
   a(8,2*j-1) = 8.0; a(8,2*j) = 1.0;
   a(9,2*j-1) = 6.0; a(9,2*j) = 2.0;
   a(10,2*j-1)= 7.0; a(10,2*j)= 3.6;
end
c(1) = 0.1; c(2) = 0.2; c(3) = 0.2; c(4) = 0.4; c(5) = 0.4;
c(6) = 0.6; c(7) = 0.3; c(8) = 0.7; c(9) = 0.5; c(10)= 0.5;
sum = 0;
for j = 1:m;
   p = 0;
   for i = 1:4
      p = p+(x(i)-a(j,i))^2;
   end
   sum = sum+1/(p+c(j));
end
y = -sum;
%--------------------------------
elseif example==43
% Shubert function
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n =2.
% 
s1 = 0; 
s2 = 0;
for i = 1:5;   
    s1 = s1+i*cos((i+1)*x(1)+i);
    s2 = s2+i*cos((i+1)*x(2)+i);
end
y = s1*s2;
%----------------------------
elseif example==44
% Sphere function 
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 30.
% 
n = 2; %(mettre ??, 5, et 10)
sum = 0;
for j = 1:n
    sum = sum+x(j)^2; 
end
y = sum;
%-----------------------------------
elseif example==45
% Sphere function 
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 30.
% 
n = 5; %(mettre ??, 5, et 10)
sum = 0;
for j = 1:n
    sum = sum+x(j)^2; 
end
y = sum;
%-----------------------------------
elseif example==46
% Sphere function 
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 30.
% 
n = 10; %(mettre ??, 5, et 10)
sum = 0;
for j = 1:n
    sum = sum+x(j)^2; 
end
y = sum;
%-----------------------------------
elseif example==47
% Sum Squares function
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 20.
% 
n = 2;
s = 0;
for j = 1:n  
    s=s+j*x(j)^2; 
end
y = s;
%-------------------------
elseif example==48
% Sum Squares function
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 20.
% 
n = 5;
sum = 0;
for j = 1:n  
    sum=sum+j*x(j)^2; 
end
y = sum;
%-------------------------
elseif example==49
% Sum Squares function
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 20.
% 
n = 10;
s = 0;
for j = 1:n  
    s=s+j*x(j)^2; 
end
y = s;
%-------------------------
elseif example==50
% Trid function
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 10.
% 
n = 6;                    % mettre 6 et 10
s1 = 0;
s2 = 0;
for j = 1:n;
    s1 = s1+(x(j)-1)^2;    
end
for j = 2:n;
    s2 = s2+x(j)*x(j-1);    
end
y = s1-s2;
%-------------------------
elseif example==51
% Trid function
% The number of variables n should be adjusted below.
% The default value of n = 10.
% 
n = 10;                    % dim 6 et 10
s1 = 0;
s2 = 0;
for j = 1:n;
    s1 = s1+(x(j)-1)^2;    
end
for j = 2:n;
    s2 = s2+x(j)*x(j-1);    
end
y = s1-s2;
%-------------------------
elseif example==52
% Zakharov function 
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 2.
% 
n = 2;                  % mettre 2, 5, 10
s1 = 0;
s2 = 0;
for j = 1:n;
    s1 = s1+x(j)^2;
    s2 = s2+0.5*j*x(j);
end
y = s1+s2^2+s2^4;
%----------------------------
elseif example==53
% Zakharov function 
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 2.
% 
n = 5;                  % mettre 2, 5, 10
s1 = 0;
s2 = 0;
for j = 1:n
    s1 = s1+x(j)^2;
    s2 = s2+0.5*j*x(j);
end
y = s1+s2^2+s2^4;
%----------------------------
elseif example==54
% Zakharov function 
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 2.
% 
n = 10;                  % mettre 2, 5, 10
s1 = 0;
s2 = 0;
for j = 1:n;
    s1 = s1+x(j)^2;
    s2 = s2+0.5*j*x(j);
end
y = s1+s2^2+s2^4;
%----------------------------

end
return