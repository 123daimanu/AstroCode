
function [ l,b,P ] = co_trans_astro( alpha,delta,p,alphap,deltap,lp_1 )
% CO_TRANS_ASTRO transforms the coordinate from one system to other.
% alpha RA of galaxy in former coordinate of galaxy
% delta DEC of galaxy in former coordinate of galaxy
% p is the position angle in former coordinate system
% alphap is the RA of the north pole of the latter coordinate system where the c
oordinate should be transformed
% deltap is the DEC of the north pole of the latter coordinate system
% lp_1 is the RA of the north pole of former coordinate system in terms of latte
r coordinate system
% Changing the degrees to radian
con    = pi/180;
alpha  = con*alpha;
delta  = con*delta;
p      = con*p;
alphap = alphap*con;
deltap = deltap*con;
lp1    = lp*con;
lp     = lp1 - pi/2
% Defining three eular matrices
R_z1 = [cos(-lp),sin(-lp),0;
       -sin(-lp),cos(-lp),0;
        0,0,1];
R_x2 = [1,0,0;
    0,cos(pi/2-deltap),sin(pi/2-deltap);
    0,-sin(pi/2-deltap),cos(pi/2-deltap)];
R_z2 = [cos(pi/2+alphap),sin(pi/2+alphap),0;
       -sin(pi/2+alphap),cos(pi/2+alphap),0;
        0,0,1];
%Transformation matrix
A_1 = R_z1*R_x2*R_z2;
%Transformed coordinates
A   = A_1*[cos(delta)*cos(alpha);cos(delta)*sin(alpha);sin(delta)];
% A = [cos(b)*cos(l);cos(b)*sin(l);sin(b)]
%value of latitude
b = asin(A(3));
%finding longitude
l = atan2(A(2),A(1));
if l < 0
    l = l + 2*pi;
end
%-----------------------------------------------------------
%to change position angle
pos  = inv(A_1);
cosk = ([-sin(delta)*cos(alpha),-sin(delta)*sin(alpha),cos(delta)]*pos)*...
    [-sin(b)*cos(l);-sin(b)*sin(l);cos(b)];
sink = ([-sin(delta)*cos(alpha),-sin(delta)*sin(alpha),cos(delta)]*pos)*...
    [-sin(l);cos(l);0];
k    = atan2(sink,cosk);
if k < 0
    k = k + 2*pi;
end
P = k + p;
if P > 2*pi
    P = P - 2*pi;
end
%--------------------------------------------------------------------
%Changing to degrees
l = l/con;
b = b/con;
P = P/con;
end
