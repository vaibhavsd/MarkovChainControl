function dxdt = func( tm,x,s,t,states,inputs,K,Xeq,B,genmat )

uelem= K*(x-Xeq);
mat= zeros(states);

for i=1:inputs
   mat= mat+ uelem(i)*genmat(:,:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proportional Controller
%     edges= [s;t]';
%     Pr= zeros(states);
%     for cnti=1:states
%         for cntj= 1:states
%             edge= [cnti,cntj];
%             if (sum(ismember(edges,edge, 'rows')))
%                 Pr(cnti,cntj)= Xeq(cnti)*x(cntj)- Xeq(cntj)*x(cnti);
%             else Pr(cnti,cntj)= 0;
%             end
%             
%         end
%     end
%     Pr(1:states+1:end)=0;
%     for newcnt=1:states
%         Pr(newcnt,newcnt)=-sum(Pr(:,newcnt));
%     end
%     mat=Pr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open Loop Control
% L= [1 -1 0 0; -1 2 -1 0; 0 -1 2 -1; 0 0 -1 1];
% U= diag([10 10 10 1/0.7])/10;
% mat= -L*U;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dxdt= mat*x;
% plot(x)
% pause(0.1);
% dxdt=B*K*x;
end


% plot(1:25,dxdt,1:25,x)
% plot(x)