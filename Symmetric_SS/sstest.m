%Paramaters 

rho = 1.225; 
theta_0 = 0;
m = 4157.174; %kg 
W = m*9.81; 
V = 300; 

%dimensions 
S = 30; %m^2
c = 2.0569; %m
b = 15.911; %m 

%aerodynamic coefficients 
C_D_0     = 0.04; 
C_L_alpha = 5.084;
e         = 0.8; 
u_c       = m/(rho*S*c);

%inertia 
K_xx = (0.019)^0.5;
K_yy = (1.3925)^0.5;
K_zz = (0.042)^0.5; 
K_xz = 0.002; 

%long. force deriv. 
C_x_u       = -0.0279;
C_x_alpha   = -0.4797;
C_x_alphad  =  0.0833;
C_x_q       = -0.2817;
C_x_delta_e = -0.0373;
C_x_0       = 2*W*sin(theta_0)/(rho*V^2*S);

%normal force deriv. 
C_z_u       = -0.3762;
C_z_alpha   = -5.7434;
C_z_alphad  = -0.0035;
C_z_q       = -5.6629;
C_z_delta_e = -0.6961;
C_z_0       = -W*2*cos(theta_0)/(rho*V^2*S);

%pitch moment deriv. 
C_m_0       =  0.0297;
C_m_u       =  0.0699;
C_m_alpha   = -0.5626;
C_m_alphad  =  0.1780;
C_m_q       = -8.7941;
C_m_delta_e = -1.1642;
C_m_T_c     = -0.0064; 

%Matrix inputs 
x_u       = V/c * C_x_u/(2*u_c);
x_alpha   = V/c * C_x_alpha/(2*u_c);
x_theta   = V/c * C_z_0/(2*u_c);

z_u       = V/c * C_z_u/(2*u_c-C_z_alphad);
z_alpha   = V/c * C_z_alpha/(2*u_c-C_z_alphad);
z_theta   = V/c * C_x_0/(2*u_c-C_z_alphad);
z_q       = V/c * (2*u_c+C_z_q)/(2*u_c-C_z_alphad);

m_u       = V/c * (C_m_u + C_z_u*(C_m_alphad)/(2*u_c-C_z_alphad))/(2*u_c*(K_yy)^2); 
m_alpha   = V/c * (C_m_alpha + C_z_alpha*(C_m_alphad)/(2*u_c-C_z_alphad))/(2*u_c*(K_yy)^2);
m_theta   = V/c * (C_x_0*(C_m_alphad)/(2*u_c-C_z_alphad))/(2*u_c*(K_yy)^2);
m_q       = V/c * (C_m_q + C_m_alphad*(2*u_c+C_z_q)/(2*u_c-C_z_alphad))/(2*u_c*(K_yy)^2);


x_delta_e = V/c * C_x_delta_e/(2*u_c);
z_delta_e = V/c * C_z_delta_e/(2*u_c-C_z_alphad);
m_delta_e = V/c * (C_m_delta_e + C_z_delta_e*(C_m_alphad)/(2*u_c-C_z_alphad))/(2*u_c*(K_yy)^2);




%State Space Matrices 

P = c/V*[-2*u_c 0                   0 0           
          0     C_z_alphad-2*u_c    0 0
          0     0                  -1 0
          0     C_m_alphad          0 -2*u_c*K_yy^2];

Q = [-C_x_u  -C_x_alpha -C_z_0  0
     -C_z_u  -C_z_alpha  C_x_0  -(C_z_q+2*u_c)
      0        0         0      -1
     -C_m_u   -C_m_alpha 0      -C_m_q];

R = [-C_x_delta_e
     -C_z_delta_e
      0
     -C_m_delta_e];
 
 A = inv(P)*Q;

 B = inv(P)*R;
 
F = [x_u x_alpha x_theta 0
     z_u z_alpha z_theta z_q
     0   0       0       V/c
     m_u m_alpha m_theta m_q];
 
E = [x_delta_e
     z_delta_e
     m_delta_e];
 
 C = [1 0 0 0
      0 1 0 0 
      0 0 1 0
      0 0 0 1];
  
 D = [0
      0
      0
      0];
 
sys = ss(A,B,C,D);

t = 0:0.01:3;
y = impulse(sys,t);
clf()
plot(t, y(:,2))
hold on
%plot(t, y(:,1))
