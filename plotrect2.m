function plotrect2(A,B,a,b)
global g
n=length(A);
if n==1
x=[A B];
y=[0 0];
plot(x,y,'b')
hold on
plot(x,y,'.w')
axis off

elseif n==2 
    hold on
x=[A(1) B(1) B(1) A(1) A(1)];
y=[A(2) A(2) B(2) B(2) A(2)];
plot(x,y,'-c')

x=[a(1) b(1)];
y=[a(2) b(2)];
%plot(x,y,':b')


plot(a(1),a(2),'.g','MarkerSize',5)
plot(b(1),b(2),'.r','MarkerSize',5)
aa=B(1)-A(1);
if a(1)==A(1)
p=a;
text(p(1)+.08*aa,p(2),['\color{green}' num2str(f(g(p)),'%4.2f')],'FontSize',6)
p=b;
text(p(1)-.08*aa,p(2),['\color{red}' num2str(f(g(p)),'%4.2f')],'FontSize',6,'horizontalAlignment','right')
else
p=a;
text(p(1)-.08*aa,p(2),['\color{green}' num2str(f(g(p)),'%4.2f')],'FontSize',6,'horizontalAlignment','right')
p=b;
text(p(1)+.08*aa,p(2),['\color{red}' num2str(f(g(p)),'%4.2f')],'FontSize',6)
end    
axis square
%xlabel('x_1')
%ylabel('x_2')
end


