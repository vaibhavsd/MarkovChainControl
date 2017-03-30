clear; clc; close all;

sys= systems(7);

dt= sys.dt;
itr= sys.itr;
particles= sys.particles;
states= sys.states;
shape= sys.shape;
s= sys.s;
t= sys.t;
Xeq= sys.Xeq;
names= sys.names;
nodes= sys.nodes;
conc= sys.conc;

n= length(s)
weights = ones(1,n)
G= digraph(s,t,weights)
figure; hold on;
hndl= plot(G,'Layout','force','EdgeLabel',G.Edges.Weight)
rotate(hndl,[0 0 1],45)
title('The Nine Node Graph Network');
hold off;
close gcf;


genmat= zeros(states,states,n);
Beq= zeros(states,n);

for i=1:n
    genmat(s(i),s(i),i)=-1;
    genmat(t(i),s(i),i)=1;
    Beq(:,i)= genmat(:,s(i),i);
end

eqfactor= zeros(1,n);

for i=1:n
    eqfactor(1,i)= Xeq(s(i));
end
eqfactor= repmat(eqfactor,states,1);
B= Beq.*eqfactor;

eta=0.1;
[T,S,V]= svd(Beq);
Ahat= zeros(states);
Ahat(states,states)= -eta;
A= T*Ahat*inv(T)
C= eye(states);
D= zeros(states,n);

K= controllerdesign(A,B,C,D,Beq,n,states)
eig(A+B*K)
eig(B*K)

tspan=[0:dt:itr*dt];
Xinit=conc;
% Xinit = (conc+1)
% Xinit = Xinit./(sum(Xinit))

[tm,y]= ode45(@(tm,y) func(tm,y,s,t,states,n,K,Xeq,B,genmat), tspan, Xinit);
figure; hold on; plot(y); title('Non-Linear System')

%y= y*particles;
% [distr, route] = particlesime3(itr,particles,nodes,conc,s,t,Xeq,K,states,dt,B,genmat);
u= K*(y'-Xeq);
figure; hold on; plot(u','LineWidth',2); title('Inputs')

% t= 0:1000;
% y2= B*u*t;
% figure; hold on; plot(y2'); title('Linearized Sys')


distr = particlesim(itr,particles,nodes,conc,Xeq,K,dt,B,genmat);


figure; hold on;

xlabel('Unit Time','Interpreter','LaTex','FontSize',14,'FontWeight','bold')
ylabel('Fraction of Agents','FontSize',14,'FontWeight','bold','Interpreter','LaTex');
% title('STOCHASTIC SIMULATION','FontSize',14,'FontWeight','bold');

linecolor={ [0    0.4470    0.7410],...
            [0.8500    0.3250    0.0980],...
            [0.9290    0.6940    0.1250],...
            [0.4940    0.1840    0.5560],...
            [0.4660    0.6740    0.1880],...
            [0.3010    0.7450    0.9330],...
            [0.6350    0.0780    0.1840] };
mx= states.*(states<7) + 6.*(states>=7);
labels= cell(1,mx);
for cnt=1:mx
    plot(tspan,y(:,cnt),'Color',linecolor{cnt},'LineWidth',1.5)
    hnd(cnt)= plot(tspan,distr(:,cnt)/particles,'Color',linecolor{cnt});
    labels{cnt}= ['x_{', num2str(cnt),'}'];
end

%hnd(cnt+1)=plot(tspan,distr(:,states+1),'--','LineWidth',1.5);
%labels{cnt+1}= ['Total Agents'];
legend(hnd,labels,'Location','east','FontSize',14)
ylim([0 0.8]);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grids

x = 0:1:shape(1);
y = 0:1:shape(2);

[X,Y] = meshgrid(x,y);
Z=0*X;
C= 0*X;
hfig= figure(1);hold on;
set(hfig, 'Position', [200 50 1000 1400]);
surf(X,Y,Z);
% xlabel('X','FontSize',16,'FontWeight','bold','Interpreter','LaTex');
% ylabel('Y','FontSize',16,'FontWeight','bold','Interpreter','LaTex');
%zlabel('Z');
% plot([0,0],[0,3],'k', 'LineWidth', 2)
% plot([1,1],[0,3],'k', 'LineWidth', 2)
% plot([2,2],[0,3],'k', 'linewidth', 2)
% plot([3,3],[0,3],'k', 'LineWidth', 2)
% 
% plot([0,3],[0,0],'k', 'linewidth', 2)
% plot([0,3],[1,1],'k', 'linewidth', 2)
% plot([0,3],[2,2],'k', 'linewidth', 2)
% plot([0,3],[3,3],'k', 'linewidth', 2)
set(gca,'xtick',[])
set(gca,'ytick',[])
pause(0.05);
for i=1:itr

listx= repmat(0:shape(1)-1,shape(2),1)';
listx= listx(:);
listy= repmat(0:shape(2)-1,shape(1),1);
listy= listy(:);

x= zeros(states,1);
for counter=1:states
    x(counter)= distr(i,counter);
    ax= listx(counter);
    bx=ax+1;
    ay= listy(counter);
    by=ay+1;
    
pxval = ax + (bx-ax)*rand(x(counter),1);
pyval = ay + (by-ay)*rand(x(counter),1);
x1hat= [1;cumsum(x)+1];
x2hat= cumsum(x);
stack(x1hat(counter):x2hat(counter),1:2)= [pxval, pyval]; 

end
pause(0.01);
if exist('hndl','var')
   delete(hndl)
end

hndl= scatter(stack(:,1),stack(:,2),'k','o','filled','LineWidth',0.01);
% txt1 = ['Itr No- ', num2str(i)];
% text(2,2.5,txt1,'Color','red','FontSize',16,'FontWeight','bold')
if(i==1)
    pause(0.5);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plotting the route of robots
% % line2= route(44,:);
% figure; hold on;
% t=1:1001;
% % plot(t,line1(1:180),'LineWidth',1.5)
% plot(t,route(21,:),'LineWidth',1.5)
% xlabel('Unit Time','FontSize',14,'FontWeight','bold','Interpreter','LaTex');
% ylabel('States','FontSize',14,'FontWeight','bold','Interpreter','LaTex')
% legend('Open Loop','Location','north')

% lmi2= route(44,:);
% figure; hold on;
% % plot(1:180,lmi1(1:180),'r','LineWidth',1.5)
% plot(1:180,ol2(1:180),'r','LineWidth',1.5)
% % plot(1:180,ol2(1:180),'g','LineWidth',1.5)
% plot(1:180,p1(1:180),'g','LineWidth',1.5)
% % plot(1:180,p1(1:180),'b','LineWidth',1.5)
% plot(1:180,lmi2(1:180),'b','LineWidth',1.5)
% xlabel('Unit Time','FontSize',14,'FontWeight','bold','Interpreter','LaTex');
% ylabel('States','FontSize',14,'FontWeight','bold','Interpreter','LaTex')
% legend('Case-1','Case-2','Case-3','Location','southeast')