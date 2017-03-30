function distr= particlesim(itr,particles,nodes,conc,Xeq,K,dt,B,genmat)

states= length(conc);
x= transpose(conc);
n= length(K(:,1));
distr= zeros(states, itr);
distr(:,1)= round(particles*x);

for i=1:itr

    
    % What did I do here?
    Xm= repmat(x,1,states);
    eps= 0.0001;
    Xn= transpose(Xm) + eps;
    
    elem= K*(x-Xeq);
    
    mat= zeros(states);
    for cnt=1:n
        mat= mat+ elem(cnt)*genmat(:,:,cnt);
    end
    

    fmat= (mat>0).*mat - (mat'<0).*(mat').*Xm./Xn;
    fmat(logical(eye(size(fmat)))) = 0;
    fmat(1:states+1:end) = -sum(fmat);
    
    fmat= (mat>0).*mat - (mat'<0).*(mat').*Xm./Xn;
    fmat(logical(eye(size(fmat)))) = 0;
    fmat(1:states+1:end) = -sum(fmat);
    
    transitionmat= fmat;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % What if controller gives unreasonably high values??
    
    maxvalue= 1/dt;
    colmatn= transitionmat;
    colmatn(logical(eye(size(fmat)))) = 0;
    sumcoln= sum(colmatn);
    for cnt=1:states
        if sumcoln(cnt)>maxvalue
            transitionmat(:,cnt)= transitionmat(:,cnt)/ sumcoln(cnt);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    transitionmat= transitionmat*dt;
    transitionmat= transitionmat';
    
    for j=1:states
        transitionmat(j,j)= 1+ transitionmat(j,j);
    end
    probmat= transitionmat;
    
    
    alpha= zeros(states,states);
    for cnt=1:states
        edges= [0 cumsum(probmat(cnt,:))];
        N= rand(round(x(cnt)*particles),1);
        alpha(:,cnt)= histcounts(N,edges)'; 
    end
    x= sum(alpha,2)/particles;
     
    distr(:,i+1)=round(x*particles);

end

distr= distr';
disp('Initial distribution: '); disp(particles*conc);
disp('Final distribution is: '); disp(distr(i+1,:));

end