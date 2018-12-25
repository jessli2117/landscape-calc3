%Project #2
%creating the equation
a1=10;
a2=5;
b1=12;
b2=9;
p=4*(a1+a2+b1+b2);
syms x y;
fx1=p/(((3*x+a1)^2)+((3*y+b1)^2)+(a1+b1));
fx2=(2*p)/(((3*x-a2)^2)+((3*y+b2)^2)+(a2+b2));
fx3=(3*p)/(((3*x+a1)^2)+((3*y-b2)^2)+(a1+b2));
fx4=(4*p)/(((3*x-a2)^2)+((3*y-b1)^2)+(a2+b1));
f=fx1-fx2-fx3+fx4;

%graphs the function
fsurf(f,[-10 10 -10 10]);
title('Land Layout');
xlabel('X-Dim');
ylabel('Y-Dim');
zlabel('Z-Dim (Height)');

%First Derivative and critical points
fx = diff(f,x);
fy = diff(f,y);
assume(x, 'real');
assume(y, 'real');
sol = vpasolve([fx, fy], [x, y]);
extremas = vpa([sol.x,sol.y],6);
xSol = sol.x;
ySol = sol.y;
%Second Derivative - final check of local maximums and minimums
fxx = diff(f,x,x);
fxy = diff(f,x,y);
fyy = diff(f,y,y);
%Make D
d = (fxx*fyy)-fxy.^2;
%Find maximums and minimums
lmax1=[];
lmax2=[];
lmin1=[];
lmin2=[];
lsaddle1=[];
lsaddle2=[];
lsaddle3=[];
for inC = 1:numel(xSol)
    numx = xSol(inC);
    numy = ySol(inC);
    newD = subs(d,[x,y],[numx,numy]); %checking d
    newFXX = subs(fxx,[x,y],[numx,numy]);%checking max/min
    if(newD>0) %if d is positive
        if(newFXX<0) %if fxx is negative (max)
            if(isempty(lmax1))
                lmax1 = [numx,numy];
            else
                lmax2 = [numx,numy];
            end
        elseif(newFXX>0) %if fxx is position (min)
            if(isempty(lmin1))
                lmin1 = [numx,numy];
            else
                lmin2 = [numx,numy];
            end
        end
    elseif(newD<0) %the saddle points
        if(isempty(lsaddle1))
            lsaddle1 = [numx,numy];
        elseif(isempty(lsaddle2))
            lsaddle2 = [numx,numy];
        else
            lsaddle3 = [numx,numy];
        end
    end
end

%Absolute Extrema
%get boundary coordinates
fxten=vpasolve(subs(f,y,10));
fxten=fxten(1);
fxnegten=real(vpasolve(subs(f,y,-10)));
fxnegten=fxnegten(1);
fyten=real(vpasolve(subs(f,x,10)));
fyten=fyten(1);
fynegten=real(vpasolve(subs(f,x,-10)));
fynegten=fynegten(1);

%try all of the coordinates
absMax =[];
absMax2 =[];
absMin =[];
absMin2 = [];
b1 = subs(f,[x,y],[fxten, 10]);
b2 = subs(f,[x,y],[fxnegten,-10]);
b3 = subs(f,[x,y],[10,fyten]);
b4 = subs(f,[x,y],[-10,fxnegten]);
%boundary max
if(b1>b2 && b1>b3 && b1>4)
    absMax = [fxten, 10,b1];
elseif(b2>b1 && b2>b3 && b2>b4)
    absMax = [fxnegten,-10,b2];
elseif(b3>b1 && b3>b2 && b3>b4)
    absMax = [10,fyten,b3];
elseif(b4>b1 && b4>b2 && b4>b3)
    absMax = [-10,fxnegten,b4];
end

%boundary min
if(b1<b2 && b1<b3 && b1<4)
    absMin = [fxten, 10,b1];
elseif(b2<b1 && b2<b3 && b2<b4)
    absMin = [fxnegten,-10,b2];
elseif(b3<b1 && b3<b2 && b3<b4)
    absMin = [10,fyten,b3];
elseif(b4<b1 && b4<b2 && b4<b3)
    absMin = [-10,fxnegten,b4];
end


%all local maximums
lmaxz1 = subs(f,[x,y],lmax1);
lmaxz2 = subs(f,[x,y],lmax2);
%find the smallest local min, then check with boundary min
if(lmaxz1>absMin(3))
    absMax = [lmax1,lmaxz1];
    if(lmaxz1<lmaxz2)
        absMax = [lmax2,lmaxz2];
        absMax2 = [lmax1,lmaxz1];
    end
end
%all local minimums
lminz1 = subs(f,[x,y],lmin1);
lminz2 = subs(f,[x,y],lmin2);
%find the smallest local min, then check with boundary min
if(lminz1<absMin(3))
    absMin = [lmin1,lminz1];
    if(lminz1>lminz2)
        absMin = [lmin2,lminz2];
        absMin2 = [lmin1,lminz1];
    end
end
%z-value of saddle points
lsaddle1z = subs(f,[x,y],lsaddle1);
saddleCoor1 = [lsaddle1, lsaddle1z];
lsaddle2z = subs(f,[x,y],lsaddle2);
saddleCoor2 = [lsaddle2, lsaddle2z];

    
%compare the maxs and minimums
disp('Absolute Mins');
disp(vpa(absMin,3));
disp('Absolute Maxs');
disp(absMax,3);

%Gradient + Directional Derivative
%Levels
fcontour(f,[-10 10 -10 10]);
%Gradient
grad = gradient(f);
grad = subs(grad,[x,y],{absMax(1),absMax(2)});
disp('gradient');
disp(grad);
%Directional Derivative
vectX = [1 0];
directX = dot(grad,vectX);
disp('directional dx');
disp(directX);
vextY = [0 1];
disp('directional dy');
directY = dot(grad,vectX);
disp(directY);

%Average Elevation
z1=@(x,y) rdivide(144,(((x.*3+10).^2)+((y.*3+12).^2)+(22)));
z2=@(x,y) rdivide(288,(((x.*3-5).^2)+((y.*3+9).^2)+(14)));
z3=@(x,y) rdivide(432,(((x.*3+10).^2)+((y.*3-9).^2)+(19)));
z4=@(x,y) rdivide(576,(((x.*3-5).^2)+((y.*3-12).^2)+(17)));
din1 = integral2(z1,-10,10,-10,10);
din2 = integral2(z2,-10,10,-10,10);
din3 = integral2(z3,-10,10,-10,10);
din4 = integral2(z4,-10,10,-10,10);
dinT = din1-din2-din3+din4;
disp(avgEle);
 