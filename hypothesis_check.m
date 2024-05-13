function [normality,homoschedasticity] = hypothesis_check(data_phase1, data_phase2, data_phase3)
    %The inputs are the data on which the statistical analysis is performed, divided in the three phases of the Stroop
    %Color Word Test. 
    %The function evaluates whether the hypothesis for the ANOVA test are
    %respected (0 no, 1 yes).
    
    normality = kstest([data_phase1(:)' data_phase2(:)' data_phase3(:)']);
    
    var_1 = var(data_phase1(:));
    var_2 = var(data_phase2(:));
    var_3 = var(data_phase3(:));
    var_array = [var_1 var_2 var_3];
    
    if all([var_1==var_2 var_2==var_3, var_1==var_3])
        homoschedasticity = 1;
    elseif max(var_array) < 2*min(var_array)
        homoschedasticity = 1;
    else
        homoschedasticity = 0;
    end
    
    %normality = 0;
end