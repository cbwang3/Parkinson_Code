classdef Trainer < handle
    
    properties (Access = private)
        classifier;
    end
    
    methods (Access = public)
        function obj = Trainer()
        end

        function train(obj,table)
            predictors = table(:,1:end-1);
            response = table(:,end);
            template = templateSVM(...
                'KernelFunction', 'linear', ...
                'KernelScale', 'auto', ...
                'Standardize', true);
            
            obj.classifier = fitcecoc(...
                predictors, ...
                response, ...
                'Learners', template, ...
                'Coding', 'onevsone');
        end
        
        function labels = test(obj,table)
            predictors = table(:,1:end-1);
            labels = predict(obj.classifier,predictors);
        end
    end
end