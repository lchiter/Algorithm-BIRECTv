function plotrect(A,B,a,b)
n=length(A);
if n==1
x=[A B];
y=[0 0];
plot(x,y,'+-b')
hold on
plot(a,0,'.g')
plot(b,0,'.r')
axis off
elseif n==2   
x=[A(1) B(1) B(1) A(1) A(1)];
y=[A(2) A(2) B(2) B(2) A(2)];
plot(x,y,'-c')
hold on

x=[a(1) b(1)];
y=[a(2) b(2)];
plot(x,y,':b')
 p=a;
plot(p(1),p(2),'.g','MarkerSize',5)
p=b;
plot(p(1),p(2),'.r','MarkerSize',5)
axis square
%xlabel('x_1')
%ylabel('x_2')
elseif n==3
x=[A(1) B(1) B(1) A(1) A(1) A(1) B(1) B(1) B(1) B(1) B(1) B(1) A(1) A(1) A(1) A(1)];
y=[A(2) A(2) B(2) B(2) A(2) A(2) A(2) A(2) A(2) B(2) B(2) B(2) B(2) B(2) B(2) A(2)];
z=[A(3) A(3) A(3) A(3) A(3) B(3) B(3) A(3) B(3) B(3) A(3) B(3) B(3) A(3) B(3) B(3)];
plot3(x,y,z,'r')
hold on
x=[a(1) b(1)];
y=[a(2) b(2)];
z=[a(3) b(3)];
plot3(x,y,z,':b')
p=a+1/3*(b-a);
plot3(p(1),p(2),p(3),'*r')
p=a+2/3*(b-a);
plot3(p(1),p(2),p(3),'*g')
axis square
xlabel('x_1')
ylabel('x_2')
zlabel('x_3')
end


