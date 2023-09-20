function [A1 B1 a1 b1 L1 U1]=Division(A,B,af,bf,lf,uf)
global  g data_ex data_fex f_eval cont
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% function [A1 B1 a1 b1 L1 U1]=Division(A,B,af,bf,lf,uf)
% global  g data_ex data_fex f_eval cont
%  Authors : Naziheddine Belkacem    (naziheddine@yahoo.fr)
% and : Lakhdar Chiter   (lchiter@univ-setif.dz)
% Created : 0510/2021
% Purpose    : Partitioning strategy 
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%
  [c i]=max(abs(B-A));
   %%%
  BB1=B;BB1(i)=BB1(i)-c/2 ;
  AA2=A;AA2(i)=AA2(i)+c/2 ;
  A1=[A AA2];
  B1=[BB1 B];
%%%%%%%%%%%  %%%%%%%%%%%%%%%%   %%%%%%%%%%%%  %%%%%%%%
af1=af;bf1=bf;
af2=af;bf2=bf;
      if bf(i)>af(i)
       af2(i)=af(i)+c/3;
       bf1(i)=bf(i)-c;%ex
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       dat=0;
       for k=1:size(data_ex,2)
           if norm(data_ex(:,k)-bf1)<1e-5
               fgbf1=data_fex(k);
               dat=1;
               f_eval=f_eval+1;
               cont=cont+1;
               break;    
           end
       end
       if dat==0
           data_ex=[data_ex bf1];
           fgbf1= f(g(bf1));
           data_fex=[data_fex fgbf1];
           f_eval=f_eval+2;
           cont=cont+2;
       end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
       L1=[lf f(g(af2))];U1=[fgbf1 uf];
      else
       af1(i)=af(i)-c/3;
       bf2(i)=bf(i)+c;%ex
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       dat=0;
       for k=1:size(data_ex,2)
           if norm(data_ex(:,k)-bf2)<1e-5
               fgbf2=data_fex(k);
               dat=1;
               f_eval=f_eval+1;
               cont=cont+1;
               break;    
           end
       end
       if dat==0
           data_ex=[data_ex bf2];
           fgbf2= f(g(bf2));
           data_fex=[data_fex fgbf2];
           f_eval=f_eval+2;
           cont=cont+2;
       end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
       L1=[f(g(af1)) lf];U1=[uf fgbf2];
      end
      a1=[af1 af2];
      b1=[bf1 bf2];     
end