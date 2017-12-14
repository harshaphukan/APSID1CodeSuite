function TabulateMPrime(F)
% Functin to tabulate SS and mprime information given appropriate array
f = figure('Position', [100 100 752 350]);
 t = uitable('Parent', f, 'Position', [75 75 700 300]);
 set(t,'Data',F);
 set(t,'ColumnName', {'SSnumber_Grain1','SchmidFactor_Grain1','SSnumber_Grain2','SchmidFactor_Grain2','Avg_SF','m_prime','abs(mp)'});
 


end