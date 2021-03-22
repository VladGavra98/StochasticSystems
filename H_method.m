clear

V     = 59.9;
g0    = 9.80665;
m     = 6035;
muc   = 113;
twmuc = 2*muc;
KY2   = 0.893;
c     = 2.022;
S     = 24.2;
lh    = 5.5;

%% TURBULENCE PARAMETERS

sigma_u = 2;
sigma_w = 3;
 Lg      = 150;
sigma_u = sigma_u;
sigma_w = sigma_w;

sigma_ug   = sigma_u/V;
sigma_ag   = sigma_w/V;

%% Derivatives - start(2) configuration
CX0  = 0.0000;
CZ0  =-1.250;
Cm0  =  0.0000;

CXu  =-0.2510;
CZu  =-2.5000;
Cmu  = 0.0000;

CXa  = 0.5120;
CZa  = -5.160;
Cma  = -0.430;

CXq  = 0.0000;
CZq  =-3.8600;
Cmq  = -7.040;

CXd  = 0.0000;
CZd  =-0.6238;
Cmd  = -1.553;

CXfa = 0.0000;
CZfa =-1.4700;
Cmfa = -3.750;

CZfug = 0.000;
Cmfug = -Cm0*lh/c;

CZfag= CZfa-CZq;
Cmfag=  Cmfa-Cmq;


xu   = (V/c)*(CXu/twmuc);
xa   = (V/c)*(CXa/twmuc);
xt   = (V/c)*(CZ0/twmuc);
xq   = 0;
xd   = (V/c)*(CXd/twmuc);
xug  = xu;
xfug = 0;
xag  = xa;
xfag = 0;

zu   = (V/c)*( CZu/(twmuc-CZfa));
za   = (V/c)*( CZa/(twmuc-CZfa));
zt   = (V/c)*(-CX0/(twmuc-CZfa));
zq   = (V/c)*((CZq+twmuc)/(twmuc-CZfa));
zd   = (V/c)*( CZd/(twmuc-CZfa));
zug  = zu;
zfug = (V/c)*( CZfug/(twmuc-CZfa));
zag  = za;
zfag = (V/c)*( CZfag/(twmuc-CZfa));

mu   = (V/c)*(( Cmu+CZu*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
ma   = (V/c)*(( Cma+CZa*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
mt   = (V/c)*((-CX0*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
mq   = (V/c)*(Cmq+Cmfa*(twmuc+CZq)/(twmuc-CZfa))/(twmuc*KY2);
md   = (V/c)*((Cmd+CZd*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
mug  = mu;
mfug = (V/c)*(Cmfug+CZfug*Cmfa/(twmuc-CZfa))/(twmuc*KY2);
mag  = ma;
mfag = (V/c)*(Cmfag+CZfag*Cmfa/(twmuc-CZfa))/(twmuc*KY2);


% Dynamics Matrix:
A_s = ([[xu, xa, xt, 0],
                [zu, za, zt, zq],
                [0, 0, 0, V/c],
                [mu, ma, mt, mq]]);




A_12 = ([[xug, xag, 0],
                 [zug-zfug*(c/Lg), zag, zfag*(c/V)],
                 [0,                0,            0],
                 [mug-mfug*(c/Lg), mag, mfag *(c/V)]]);



A_22 = ([[-V/Lg,      0,          0],
                 [0,          0,          1],
                 [0, -V^2/(Lg^2), -2*(V/Lg)]]);


% A = [A_s  A12]
%     [0    A22]
A = [[A_s, A_12],
     [zeros(3,4), A_22]];



% Input Mtrix:
B = ([[ xd,    0.0,                                       0.0                                   ],
                    [ zd,   zfug*(c/V)*sigma_ug*sqrt(2*V/Lg), zfag*(c/V)*sigma_ag*sqrt(3*V/Lg)],
                    [ 0.0 , 0.0,                                        0.0                                   ],
                    [  md,  mfug*(c/V)*sigma_ug*sqrt(2*V/Lg),  mfag*(c/V)*sigma_ag*sqrt(3*V/Lg)],
                    [ 0.0 , sigma_ug*sqrt(2*V/Lg),                      0.0                           ],
                    [ 0.0,          0.0   ,                               sigma_ag*sqrt(3*V/Lg)],
                    [  0.0,  0.0,                           (1-2*sqrt(3))*sigma_ag*power(V/Lg,1.5)]]);

% Output matrix:
C = ([[1, 0, 0, 0, 0, 0, 0],
                   [0, 1, 0, 0, 0, 0, 0],
                   [0, 0, 1, 0, 0, 0, 0],
                   [0, 0, 0, 1, 0, 0, 0],
                   [-V*zu/g0, -V*za/g0, -zt/g0, (V*V/c - V*zq)/g0, -V*(zug-zfug*(c/Lg))/g0, -V*zag/g0, -V*zfag*(c/V)/g0]]);
               
% Feed-thourgh matrix:
D      = zeros(5,3);
D(2,:) = -V*([zd, zfug*(c/V)*sigma_ug*sqrt(2*V/Lg), zfag*(c/V)*sigma_ag*sqrt(3*V/Lg)])/g0;




model_ss = ss(A,B,C,D);

%% Simulation Parameters:

%TIME AXIS INPUT VECTOR DEFINITION
dt    = 0.01;               % sec
T     = 10000;               % sec
t     = 0:dt:T-dt;  % sec - check for lickage
N     = length(t);             % number of samples

% Selected turbulence input:
    % 2 for w1  = horizontal
    % 3 for w3  = vertical

windex       = 3;    % CHANGE THIS

turb_lst = {'zero', "horizontal", 'vertical'};
disp('Turbulence mode: ');
disp(string(turb_lst(windex)));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method 0.                Lyapunov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wc = 1;


Bin = B(:,windex);
L = lyap(A, Bin * Wc * Bin' );

for i=1:5
    var_lyap(i) = L(i,i);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method 3.                 TIME VARIANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


u_noise(windex,1:N) = randn(1,N)/sqrt(dt);
yout = lsim(model_ss,u_noise,t);

for i=1:5
    y = yout(1:N,i);
    var_time(i) = var(y,0,1);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method 2. USING THE IMPULSE RESPONSE METHOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ZERO INPUT, INITIAL CONDITION EQUALS B (for input 3 in this case)
u = zeros(3,N); x0=B(:,windex);

% CALCULATION OF IMPULSE RESPONSES
h = lsim(A,B,C,D,u,t,x0);	

% PLOT IMPULSE RESPONSE
subplot(2,1,1)
plot(t,h(:,1)); xlabel('Time [sec]'); ylabel('h u w3(t)');
subplot(2,1,2)
plot(t,h(:,2)); xlabel('Time [sec]'); ylabel('h alpha w3(t)');


% CALCULATION OF PRODUCT MATRIX OF IMPULSE RESPONSES
h11=h(:,1).*h(:,1);h12=h(:,1).*h(:,2);h13=h(:,1).*h(:,3);h14=h(:,1).*h(:,4);h15=h(:,1).*h(:,5);
                   h22=h(:,2).*h(:,2);h23=h(:,2).*h(:,3);h24=h(:,2).*h(:,4);h25=h(:,2).*h(:,5);
                                      h33=h(:,3).*h(:,3);h34=h(:,3).*h(:,4);h35=h(:,3).*h(:,5);
                                                         h44=h(:,4).*h(:,4);h45=h(:,4).*h(:,5);
                                                                            h55=h(:,5).*h(:,5);
 % PLOT (CROSS) PRODUCTS OF IMPULSE RESPONSES
plot(t,h11); xlabel('Time [sec]'); ylabel('h1*h1(t)');


% INTEGRATION OF PRODUCT MATRIX OF IMPULSE RESPONSES
var11(1)=0; var12(1)=0; var13(1)=0; var14(1)=0;
            var22(1)=0; var23(1)=0; var24(1)=0;
                        var33(1)=0; var34(1)=0;
                                    var44(1)=0; var55(1)=0;
                                    
dth11 = dt*h11; dth12 = dt*h12; dth13 = dt*h13; dth14 = dt*h14;
                dth22 = dt*h22; dth23 = dt*h23; dth24 = dt*h24;
                                dth33 = dt*h33; dth34 = dt*h34;
                                                dth44 = dt*h44;
                                                dth55 = dt*h55;
                                                
 % Integral loop:                                                               
for i=1:N-1
    
    var11(i+1) = var11(i) + dth11(i);
    var12(i+1) = var12(i) + dth12(i);
    var13(i+1) = var13(i) + dth13(i);
    var14(i+1) = var14(i) + dth14(i);
    var22(i+1) = var22(i) + dth22(i);
    var23(i+1) = var23(i) + dth23(i);
    var24(i+1) = var24(i) + dth24(i);
    var33(i+1) = var33(i) + dth33(i);
    var34(i+1) = var34(i) + dth34(i);
    var44(i+1) = var44(i) + dth44(i);
    var55(i+1) = var55(i) + dth55(i);
end

fprintf("var_H = %1.10f, var_lya = %1.10f var_time = %1.10f \n",var11(end), var_lyap(1), var_time(1));
fprintf("var_H = %1.10f, var_lya = %1.10f var_time = %1.10f \n",var22(end), var_lyap(2), var_time(2));
fprintf("var_H = %1.10f, var_lya = %1.10f var_time = %1.10f \n",var33(end), var_lyap(3), var_time(3));
fprintf("var_H = %1.10f, var_lya = %1.10f var_time = %1.10f \n",var44(end), var_lyap(4), var_time(4));
fprintf("var_H = %1.10f, var_lya = %1.10f var_time = %1.10f \n",var55(end), var_lyap(5), var_time(5));
