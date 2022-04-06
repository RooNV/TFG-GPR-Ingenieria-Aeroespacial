function gprMdl = Mat_gpr(Xtrain, ftrain, Kfn)

    rng(0,'twister'); % For reproducibility
    
    
    H0 = [1 0.4 85];
    sigmaN0 = 0.9178;
    gprMdl = fitrgp(Xtrain,ftrain,'KernelFunction',Kfn,'KernelParameters', ...
          H0,'Sigma',sigmaN0);
%           H0);
%       , 'Standardize', 1);
    
%     gprMdl = fitrgp(Xtrain,ftrain,'KernelFunction','squaredexponential',...
%     'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
%     struct('AcquisitionFunctionName','expected-improvement-plus'));
    
    % Parámetros
    names = gprMdl.KernelInformation.KernelParameterNames
    values_1 = gprMdl.KernelInformation.KernelParameters
    sigma = gprMdl.Sigma

end