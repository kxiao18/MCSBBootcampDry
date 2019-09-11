% horizontal axis is the ratio of kinase to phosphatase (say, in log scale), and
% the vertical axis is the fraction of activated protein in steady state.

% set up initial conditions
ITot = 1;
KTot = 1 ;
PTot = 1;
I0 = 1;
A0 = 0;
IK0=2;
AP0=1;

% set rates
kAon = 10;
kAoff = 10;
kIon = 10;
kIoff = 10;
kAcat = 100;
kIcat = 10;

% create differential equations
dAdt = @(A,I,AP,IK) -kAon*A*(PTot-AP) + kAoff*AP + kAcat*IK;
dIdt = @(A,I,AP,IK) -kIon*I*(KTot-IK) + kIoff*IK + kIcat*AP;
dAPdt = @(A,I,AP,IK) kAon*A*(PTot-AP) - kAoff*AP - kIcat*AP;
dIKdt = @(A,I,AP,IK) kIon*I*(KTot-IK) - kIoff*IK - kAcat*IK;

dxdt = @(t,x) [ dAdt(x(1),x(2),x(3),x(4));
        dIdt(x(1),x(2),x(3),x(4));
        dAPdt(x(1),x(2),x(3),x(4)) ;
        dIKdt(x(1),x(2),x(3),x(4));]
    
% solve system
[T,X] = ode45(dxdt,[0,2],[A0,I0,AP0,IK0]);

% plot
figure(1); clf; hold on; box on;
plot(T,X,'LineWidth',2);
plot(T,sum(X,2),'--k','LineWidth',2);

legend('A', 'I', 'AP','IK');