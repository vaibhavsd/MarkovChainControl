function val= init_robots(n,conc,nodes)

val= zeros(n,1);
modes= length(conc);
conc= [0, conc]
particles= n*conc;
particles= cumsum(particles);

for i=1:modes
   for j= particles(i)+1:particles(i+1)
       val(j)= nodes(i);
   end
end

end