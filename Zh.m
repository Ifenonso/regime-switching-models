function [UL,WL,YL] = Zh(K,sf,U,W,U1,W1,h,varargin)

format long
for i = 1:length(sf)
    gh(i) = K;
end
ab1 = isreal(sf); ac1 = isreal(U); ad1 = isreal(W); 
if ab1 == 0 && ac1 == 0 && ad1 == 0 
    sf = gh; U = U1; W = W1;
end

if min(sf)<0 
    sf = gh; U = U1; W = W1;
end

a = length(sf);
for i=1:a
    for j = 1:a
        b = log(sf(j)/sf(i));
        if -b==0
            Ub = K-sf(j)*exp(-b);
            Wb = -sf(j)*exp(-b);
            Yb = -sf(j)*exp(-b);
        elseif -b<0
            Ub = K-sf(j)*exp(-b);
            Wb = -sf(j)*exp(-b);
            Yb = -sf(j)*exp(-b);
%         elseif 
        else    
            ck = -b; dk = -b/h; ek = fix(dk);             
            aa = (U(j,ek+2)-U(j,ek+1))/h; ab = (aa-W(j,ek+1))/h;
            ac = (W(j,ek+2)-aa)/h; ad = (ac-ab)/h;
            ae = (U(j,ek+3)-U(j,ek+2))/h; af = (ae-W(j,ek+2))/h;
            ag = (af-ac)/(2*h); ah = (ag-ad)/(2*h); 
            ai = (W(j,ek+3)-ae)/h; aj = (ai-af)/h;
            ak = (aj-ag)/2*h; am = (ak-ah)*2*h;    
            
            Ub = U(j,ek+1)+W(j,ek+1)*(ck-ek*h)+ab*((ck-ek*h)^2)+ad*...
                ((ck-ek*h)^2)*(ck-(ek+1)*h)+ah*((ck-ek*h)^2)*((ck-...
                (ek+1)*h)^2)+am*((ck-ek*h)^2)*((ck-(ek+1)*h)^2)*...
                (ck-(ek+2)*h);            
            
            Wb = W(j,ek+1)+ab*2*(ck-ek*h)+ad*(((ck-ek*h)^2)+2*(ck-...
                ek*h)*(ck-(ek+1)*h))+ah*(2*(ck-ek*h)*((ck-(ek+1)*h)^2)...
                +2*((ck-ek*h)^2)*(ck-(ek+1)*h))+am*(2*(ck-ek*h)*((ck-...
                (ek+1)*h)^2)*(ck-(ek+2)*h)+2*((ck-ek*h)^2)*(ck-(ek+1)...
                *h)*(ck-(ek+2)*h)+((ck-ek*h)^2)*((ck-(ek+1)*h)^2)); 
            
            Yb = 2*ab+ad*(4*(ck-ek*h)+2*(ck-(ek+1)*h))+ah*(2*((ck-(ek+1)...
                *h)^2)+8*(ck-(ek+1)*h)*(ck-ek*h)+2*((ck-ek*h)^2))+am*(4*...
                ((ck-ek*h)^2)*(ck-(ek+1)*h)+4*((ck-(ek+1)*h)^2)*(ck-...
                ek*h)+2*((ck-ek*h)^2)*(ck-(ek+2)*h)+2*((ck-(ek+1)*h)^2)*...
                (ck-(ek+2)*h)+8*(ck-ek*h)*(ck-(ek+1)*h)*(ck-(ek+2)*h));
        end            
        UL((j+a*(i-1))) = Ub;
        WL((j+a*(i-1))) = Wb;
        YL((j+a*(i-1))) = Yb;
    end
end
end
