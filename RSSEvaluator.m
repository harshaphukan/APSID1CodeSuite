%Function to Evaluate Resolved shear Stress
% Inputs: Slip/Twin Plane normal(n) (Unit Vector),Shear/Slip Direction(m)
function RSS=RSSEvaluator(m,n,STens)


Sigma=[STens(1),STens(4),STens(5);
       STens(4),STens(2),STens(6);
       STens(5),STens(6),STens(3)];
% Take dyadic product of m and n
SM=0.5*(m'*n+n'*m);
RSS=trace(Sigma*SM);

end