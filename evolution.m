function route= evolution(val,probmat,nodes)
    
    n=length(val(:,1));
    route= zeros(n,2);
    route(:,1)= val;  
    pmat=(cumsum(probmat'))';
    [~, c]= size(pmat);
    pmat= [zeros(c,1), pmat];
    r= length(nodes);
    p= rand(n,1);
        
    for j=1:n %i= no of particles; 
       
        for count=1:r
            %nodes(count)
        if route(j,1)== nodes(count)
            prob= pmat(count,:);
        end
        end
        

        % disp(p(j));
        for k=1:length(prob)-1  %j= no of cols; 
            
        if p(j)>prob(k) && p(j)<=prob(k+1)
            %do something
            route(j,1+1)= nodes(k);
        end
        end
                
    end
end