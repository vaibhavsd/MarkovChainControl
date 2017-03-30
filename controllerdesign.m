function K= controllerdesign(A,B,C,D, Beq,n,states)

val=4;
eta=0.1;
Kstruct= (Beq'>0).*Beq' - (Beq'<0).*Beq';
Phat= sdpvar(states,states);
mat5= diag(ones(states,1));
P= Phat.*mat5;
Zhat= sdpvar(n,states);
Z= Zhat.*Kstruct;


switch(val)
    
    case 1
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stabilizability Controller Design 
constraints= [P>=eta*eye(states); A*P+B*Z+P*A'+Z'*B'<= -eta*eye(states);];
optimize(constraints);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    case 2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D- Stability Controller Design 
tr= 3;
r= 1.8/tr;
alpha= 1; 
c=-1;


mat1=[-r*P, A*P+B*Z; (A*P+B*Z)', -r*P];
mat2=[c*(A*P+B*Z+(A*P+B*Z)'), (A*P+B*Z-(A*P+B*Z)');...
      (A*P+B*Z-(A*P+B*Z)')', c*(A*P+B*Z+(A*P+B*Z)')];
  
c1= [P>=eta*eye(states)];
c2= [mat1<= -eta*eye(2*states)];
c3= [A*P+B*Z+(A*P+B*Z)'+2*alpha*P<= -eta*eye(states)];
c4= [mat2<= -eta*eye(2*states)];
constraints= [c1; c2; c3; c4];
optimize(constraints);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    case 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H-infinity based controller Rodriguez
% measure numbers of inputs and outputs 
% sys= ss(A,B,C,D);
% NMEAS= size(C2,1);
% NCON =size(B2,2);
% hinfsyn(sys,NMEAS, NCON)
% 
% 
% Ap= A; Bp= B; Cp=1; Dp=0;  
% A1= 0; A2= 0; A3= 0;
% B1= 0; B2= 0; B3= 0;
% C1= 0; C2= 0; C3= 0;
% D1= 1; D2= 1; D3= 1;
% 
% A= [A1 0 0 -B1*Cp; 
%      0 A2 0 0;
%      0 0 A3 B3*Cp;
%      0 0 0 Ap];
% B1= [0 0 0 0];
% B2= [0 B2 0 Bp];
% C1= [C1 0 0 -D1*Cp;
%      0 C2 0 0;
%      0 0 C3 D3*Cp];
% C2= [0 0 0 -Cp];
% D11= [D1;0;0];
% D12= [0;D2;0];
% D21= 1;
% D22= 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    case 4
% H-infinity based controller Peet

B1= 0.01*eye(states,n);
B2= B;

C1= [eye(states);zeros(n,states)];
C2= eye(states);

D11= [zeros(states,n); zeros(n,n)];
D12= [zeros(states,n); eye(n)];
D21= zeros(states,n);
D22= zeros(states,n);

eta=.001;    % degree of strict positivity   
ns=size(A,1);   % number of states
nc=size(B2,2);  % number of actuators
nd=size(B1,2);  % number of external inputs
nr=size(C1,1);  % number of regulated outputs
C2t=eye(ns); D21t=zeros(ns,nd); D22t=zeros(ns,nc);
nm=size(C2,1);  % number of sensors


% H-infinity State Feedback Controller Synthesis


% Declare the variables
gamma=sdpvar(1);               % represents the bound on the H-infinity norm of the CL system.
X=sdpvar(ns);
mat1= diag(ones(1,ns));
% X= X.*mat1;
X= 0.01*eye(ns);
Zhat=sdpvar(nc,ns,'full');
Z= Zhat.*Kstruct;

W=sdpvar(nr);

% declare constraints
MAT=[A*X+X*A'+B2*Z+Z'*B2'       B1              (C1*X+D12*Z)';
     B1'                        -gamma*eye(nd)   D11';
     C1*X+D12*Z                 D11              -gamma*eye(nr)];
F=[MAT<=0];
F=[F;X>=eta*eye(ns)];

OPTIONS = sdpsettings('solver','sedumi');

% Solve the LMI, minimizing gamma
optimize(F,gamma,OPTIONS);
gamman=value(gamma);
disp('The value of gamma is: '); disp(gamman);

% retrieve decision variables
P=value(X); Z=value(Z); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

P= value(P)
Z= value(Z);
K= Z*inv(P)

end