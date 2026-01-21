function [shape, deriv] = my_shape_fun(exisp,etasp,nnode)
format long;

%% 3 nodes elements: triangular finite elements
if(nnode == 3)
    % local coordinates of Gaussian integration point 
    s  = exisp;
    t  = etasp;
    p  = 1.0 - s - t;
    
    % shape function values of the triangular element
    shape(1) = p;
    shape(2) = s;
    shape(3) = t;
    
    % the derivatives of shape functions with respect to local coordinates
    deriv(1, 1) = -1.0;
    deriv(1, 2) = 1.0;
    deriv(1, 3) = 0.0;
    
    deriv(2, 1) = -1.0;
    deriv(2, 2) = 0.0;
    deriv(2, 3) = 1.0;
    
end

%% 4 nodes elements : quadrilateral finite elements
if(nnode == 4 )
    % local coordinates of Gaussian integration point
    s  = exisp;
    t  = etasp;
    st = s * t;
    
    % shape function values of the quadrilateral element
    shape(1)=(1.0-t-s+st)*0.25;
    shape(2)=(1.0-t+s-st)*0.25;
    shape(3)=(1.0+t+s+st)*0.25;
    shape(4)=(1.0+t-s-st)*0.25;
    
    % the derivatives of shape functions with respect to local coordinates
    deriv(1, 1) = (-1.0+t) * 0.25;
    deriv(1, 2) = (1.0-t)  * 0.25;
    deriv(1, 3) = (1.0+t)  * 0.25;
    deriv(1, 4) = (-1.0-t) * 0.25;
    %
    deriv(2, 1) = (-1.0+s) * 0.25;
    deriv(2, 2) = (-1.0-s) * 0.25;
    deriv(2, 3) = (1.0+s)  * 0.25;
    deriv(2, 4) = (1.0-s)  * 0.25;
    
end

%%  8 nodes elements
if(nnode == 8 )
    %  gaussian integration point coordinate information
    s=exisp;
    t=etasp;
    s2=2.0*s;
    t2=2.0*t;
    ss=s*s;
    tt=t*t;
    st=s*t;
    sst=s*s*t;
    stt=s*t*t;
    st2=2.0*s*t;
    
    % the function values of the shape functions at the Gaussian nodes in local coordinates
    shape(1)=(-1.0+st+ss+tt-sst-stt)*0.25;
    shape(2)=(1.0-t-ss+sst)*0.5;
    shape(3)=(-1.0-st+ss+tt-sst+stt)*0.25;
    shape(4)=(1.0+s-tt-stt)*0.5;
    shape(5)=(-1.0+st+ss+tt+sst+stt)*0.25;
    shape(6)=(1.0+t-ss-sst)*0.5;
    shape(7)=(-1.0-st+ss+tt+sst-stt)*0.25;
    shape(8)=(1.0-s-tt+stt)*0.5;
    
    % the derivative values
    deriv(1,1)=(t+s2-st2-tt)*0.25;
    deriv(1,2)=-s+st;
    deriv(1,3)=(-t+s2-st2+tt)*0.25;
    deriv(1,4)=(1.0-tt)*0.5;
    deriv(1,5)=(t+s2+st2+tt)*0.25;
    deriv(1,6)=-s-st;
    deriv(1,7)=(-t+s2+st2-tt)*0.25;
    deriv(1,8)=(-1.0+tt)*0.5;
    %
    deriv(2,1)=(s+t2-ss-st2)*0.25;
    deriv(2,2)=(-1.0+ss)*0.5;
    deriv(2,3)=(-s+t2-ss+st2)*0.25;
    deriv(2,4)=-t-st;
    deriv(2,5)=(s+t2+ss+st2)*0.25;
    deriv(2,6)=(1.0-ss)*0.5;
    deriv(2,7)=(-s+t2+ss-st2)*0.25;
    deriv(2,8)=-t+st;
end
end %endfunction



