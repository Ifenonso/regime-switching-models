function [UL,WL] = NewtHermite(K,sf,U,W,U1,W1,nx,h,x_max,varargin)

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
UL = zeros(a^2,nx+1);
WL = zeros(a^2,nx+1);
for i=1:a
    for j = 1:a
            for k = 1:nx+1
                b = log(sf(j)/sf(i));
                if (k-1)*h-b<=0
                    ck = (k-1)*h-b;
                    Ub(k) = K-sf(j)*exp(ck);
                    Wb(k) = -sf(j)*exp(ck);
                elseif (k-1)*h-b>=x_max
                    Ub(k) = 0;
                    Wb(k) = 0;
                else
                    ck = (k-1)*h-b; dk = ((k-1)*h-b)/h;
                    ek = fix(dk); 
                    
                    aa = (U(j,ek+2)-U(j,ek+1))/h; ab = (aa-W(j,ek+1))/h;
                    ac = (W(j,ek+2)-aa)/h; ad = (ac-ab)/h;
                    
                    Ub(k) = U(j,ek+1)+W(j,ek+1)*(ck-ek*h)+ab*...
                        ((ck-ek*h)^2)+ad*((ck-ek*h)^2)*(ck-(ek+1)*h);
                    Wb(k) = W(j,ek+1)+ab*2*(ck-ek*h)+ad*(((ck-ek*h)^2)+...
                        2*(ck-ek*h)*(ck-(ek+1)*h));
                end
                
            end
            WL((j+a*(i-1)),:) = Wb;
            UL((j+a*(i-1)),:) = Ub;
    end
        
end

end