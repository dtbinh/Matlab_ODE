function x=Lorenz(x0,p,opcion)


if 0>1
    clc,clear all,close all
    x0=[-8,8,27];
    tspan=[0,20];
    %Lorenz equations.
    % sig=p(1);
    % beta=p(2);
    % rho=p(3);
    p=[10,8/3,28];
    
    if 0
        %opcion 1: odes
        opcion='ode45';% o cualquier ode que tenga nomenclatura odeXX
        opcion='ode23s';%...      
    else
        %opcion 2: simulink
        opcion='SIM';
    end
    
    x=Lorenz(x0,p,opcion);
    
end

% VARIABLES DE ENTRADA
% x0: initianal values of x,y,z OBLIGATORIO
% p: parameters of sig,rho,beta OBLIGATORIO
% opcion: 'SIM' or 'ODEXX'
% tspan: limits of t


if nargin<2
    error('Faltan datos.')
elseif exist('x0')==0 || exist('p')==0
    error('Faltan datos OBLIGATORIOS.')% exigencia del enunciado
end

if isempty('opción')==1
    opcion='SIM';% default
end

if exist('tspan')==0
    tspan=[0,20];% default
end

%%
% si opcion='SIM' opera por simulink
% si opcion='ODEXX' opera por odexx, por ejemplo ode45

if strcmp(opcion,'SIM')==1
    disp('SIM method')

    % http://www.mathworks.com/matlabcentral/fileexchange/28451-solution-of-differential-equations-with-matlab-simulink-lorenz-attractor-case-study
    % SimOut = sim('model', ParameterStruct);
    
    %opts=simset('MaxStep',1e-2);
    %[t,x]=sim('lorenz_sim',tspan,opts);% metodo con ctes
    
    N=1e5;
    u(:,1)=ones(N,1)*x0(1);
    u(:,2)=ones(N,1)*x0(2);
    u(:,3)=ones(N,1)*x0(3);
    u(:,4)=ones(N,1)*p(1);
    u(:,5)=ones(N,1)*p(2);
    u(:,6)=ones(N,1)*p(3);
    t=linspace(tspan(1),tspan(2),N)';
    [t,x]=sim('lorenz_sim2',t,[],[t,u]);% metodo con inputs
    
else
    % ode45(function,domain,initial condition)
    %ej. [t,x]=ode45(@lorenzFcn,tspan,x0,[],p);
        
    disp(['Resolviendo por método ',opcion,'.']);
    
    % CALCULOS
    [t,x]=eval([opcion,'(@lorenzFcn,tspan,x0,[],p)']);
        
    
    % breve respado de eval
    % ej. x=rand(1e2);a='plot',eval([a,'(x)'])
end

figure(1),plot3(x(:,1),x(:,2),x(:,3)),grid on
titleStr=sprintf('Lorenz function\nx_0=%d, y_0=%d, z_0=%d\nsigma=%1.2g, beta=%1.2g, rho=%1.2g',...
    x0(1),x0(2),x0(3),p(1),p(2),p(3));
    
title(titleStr,'fontweight','bold','fontsize',20);
figure(2),hold on
N=3;

col=hsv(N);
legendStr = cell(N,1);  
for k=1:N
    legendStr(k) = {sprintf('x(:,%d),t',k)};   
    plot(t,x(:,k),'color',col(k,:))
end
grid on,hold off
legend(legendStr,'Location','Best','fontweight','bold');
title(titleStr,'fontweight','bold','fontsize',20);  
xlabel('time (sec)');  
ylabel('amplitude');  
return

function x =lorenzFcn(t,x,p);
%LORENZ: Computes the derivatives involved in solving the
%Lorenz equations.
sig=p(1);
beta=p(2);
rho=p(3);
% dx/dt=-sig*x + sig*y=-p(1)*(u(1)+u(2))
% dy/dt=rho*x-y-x*z=p(2)*u(1)-u(2)-u(1)*u(3)
% dz/dt=-beta*z + x*y=-p(3)*u(3)+u(1)*u(2)
x=[-sig*x(1)+sig*x(2);rho*x(1)-x(2)-x(1)*x(3); -beta*x(3)+x(1)*x(2)];
% Observe that:
% x is stored as x(1)
% y is stored as x(2)
% z is stored as x(3)

return

