%--------------------------------------------------------------------------
% Selection of potential optimal hyper-rectangles
%  Authors : Naziheddine Belkacem    (naziheddine@yahoo.fr)
% and : Lakhdar Chiter   (lchiter@univ-setif.dz)
% Created : 05/10/2021
% Purpose    : Running BIRECTv or BIRECTv-l 
%--------------------------------------------------------------------------
function [ll]=po_opt(nor,D,fmin,epsilon)
ll=[];
noru=unique(nor);%noru=sort(noru,'descend');
for j=1:length(noru)
    I=find(nor==noru(j));
    [mn ind]= min(D(I));
    ind=I(ind);
%%% All hyper-rectangles having the same minimum value of f and the same size %%%%
%%% uncomment ind2=find(mn==D(I),1,'first'); and comment ind2=find(mn==D(I)); to run BIRECT-vl %%%
%%% do the opposite if you want to run BIRECT-v %%%
ind2=find(mn==D(I));
% ind2=find(mn==D(I),1,'first');
    ind2=I(ind2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I1=find(nor<noru(j));
    I2=find(nor>noru(j));
    
    if isempty(I1)
        Kmax=-inf;
    else
        Kmax=max((D(ind)-D(I1))./(noru(j)-nor(I1)));
    end
    
    if isempty(I2)
        Kmin=inf;
    else
        Kmin=min((D(I2)-D(ind))./(nor(I2)-noru(j)));
    end
    
    if Kmin>0 && Kmax<=Kmin
        %%%%%%%%%%%
        
        if fmin~=0&& epsilon<=(fmin-D(ind))/abs(fmin)+noru(j)/abs(fmin)*Kmin
            ll=[ll ind2];
        elseif fmin==0&&D(ind)<=noru(j)*Kmin
            ll=[ll ind2];
        end
%         if (D(ind)-fmin+epsilon*abs(fmin))/noru(j)<=Kmin
%           ll=[ll ind2];  
%         end
    end
    
    
    
    %         if fmin~=0
%             if epsilon<=(fmin-D(ind))/abs(fmin)+noru(j)*Kmin/abs(fmin)
%                ll=[ll ind2];
%             end
%         else
%             if D(ind)<=noru(j)*Kmin
%                 ll=[ll ind2];
%             end
%         end
         
        
end
end        
