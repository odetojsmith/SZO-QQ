import mosek.fusion.*;
 
tic
M = Model('problem');
 
z = M.variable('z',3); % new variable
z0 = Expr.constTerm(ones(3,1)); % new constant
x = Expr.add(z,z0); % x = z + z0
 
t = M.variable('t',1); % new varaible
 
M.constraint(Expr.vstack(t,x), Domain.inQCone()); % add a quadratic cone constraint
M.constraint(z.index(0), Domain.greaterThan(0.5)); % add â‰¥ constraint 
% !!!!!!! fusion uses zero-based indexing !!!!!!!!
 
M.objective('obj', ObjectiveSense.Minimize, t); % specify objective
 
M.solve(); % aaaand solve
 
 
z.level(); % get value of variable z
t.level() % get value of variable t
toc