clear;clc;

Nx = 5;     Ny = 2;

L = 0.5;    H = 0.2;    t = 0.05;

[ x , y ] = ndgrid( linspace(0 , L , Nx+1 ) , linspace( 0 , H , Ny+1 ));
% plot( x , y , 'xb' )
line(x,y,'color','k')
axis( [ -0.1 L+0.1 -0.1 H+0.1])
grid on
