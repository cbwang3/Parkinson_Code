close all;

path = genpath('./');
addpath(path);
loadClasses;

fileNames = {'1-lukas','2-dino','3-luis', '4-mauricio','5-leo','7-noa','8-alex','9-jakob','10-michi','11-roko','12-theo','13-aldin'};
%fileNames = {'9-jakob'};

global exercisesTable;
lazyLoadExercisesTableForSegmentation(fileNames);

nFeatures = size(exercisesTable,2)-1;
predictors = exercisesTable(:,1:nFeatures);
predictors = table2array(predictors);
responses = exercisesTable.label;

%{
%% Filter Approach
%prepare data
relevantInstanceIdxs = (responses == 1);
relevantInstances = predictors(relevantInstanceIdxs,:);
irrelevantInstances = predictors(~relevantInstanceIdxs,:);
[h,p,ci,stat] = ttest2(relevantInstances,irrelevantInstances,'Vartype','unequal');

%this plots indicates the percentage of features that have a low 'p' index
%low p indexes are good, they tell the likelihood that the features come from
%independent data sources
ecdf(p);
xlabel('P value');
ylabel('CDF value');
%}

%{
%% select features based on Filter approach
[~,featureIdxSortbyP] = sort(p,2); % sort the features

maxFeatures = 200;
testMCE = zeros(1,maxFeatures);
testMCEResub = zeros(1,maxFeatures);

%ytest = classify(xtest, xtrain, ytrain, 'diaglinear');

classf = @(xtrain,ytrain,xtest,ytest) ...
    sum(ytest ~= classify(xtest,xtrain,ytrain,'quadratic'));

%cv = cvpartition(responses,'kfold',10);
cv = cvpartition(responses,'holdout',1000);
cvResub = cvpartition(responses,'resubstitution');%testing = training
excludeIdxs = [16, 73, 80, 85, 87, 93, 100, 102, 103, 104, 108, 110, 111, 112, 114, 116, 120, 123, 124, 130, 131, 132, 137, 139, 140, 143, 144, 146, 147, 149, 150, 151, 152, 155, 158, 159, 162, 166, 167, 168, 174, 175, 177, 181, 183, 184, 186, 187, 188, 191, 195, 196, 197, 198];

%cv = cvpartition(responses,'holdout',1000);
for i = 1 : maxFeatures
    fs = featureIdxSortbyP(1:i);
    fs = setdiff(fs,featureIdxSortbyP(excludeIdxs));
    try
        res = crossval(classf,predictors(:,fs),responses,'partition',cv)' ./ cv.TestSize;
        resResub = crossval(classf,predictors(:,fs),responses,'partition',cvResub)' / cvResub.TestSize;
    catch ME
        fprintf('error for %d\n',i);
        excludeIdxs = [excludeIdxs, i];
        display(ME.message);
    end
    
    testMCE(i) = res;
    testMCEResub(i) = resResub;
end

nfs = [1:maxFeatures];
plot(nfs, testMCE,'o',nfs,testMCEResub,'r^');
xlabel('Number of Features');
ylabel('MCE');
legend({'MCE on the test set' 'Resubstitution MCE'},'location','NW');
title('Simple Filter Feature Selection Method');

[minVal, minIdx] = min(testMCE);

fprintf('Min MCE error: %.4f using features < %d\n',minVal, minIdx);

minIdx = 138;
fprintf('Choosing minIdx = %d with MCE error: %.4f\n',minIdx, testMCE(minIdx));


preSelectedFeatureIdxs = featureIdxSortbyP(1:minIdx);
preSelectedFeatureIdxs = setdiff(preSelectedFeatureIdxs,featureIdxSortbyP(excludeIdxs));
res = crossval(classf,predictors(:,preSelectedFeatureIdxs),responses,'partition',cv)' ./ cv.TestSize;

fprintf('Retesting error: %.4f\n',res);

fprintf('Features after quadratic filter: \n');
for i = 1 : length(preSelectedFeatureIdxs)
    feature = preSelectedFeatureIdxs(i);
    fprintf('%d, ',feature);
end
%}

%{
%% Apply cross correlation to eliminate highly correlated features
preSelectedFeatureIdxs = [1, 2, 3, 4, 5, 6, 9, 10, 11, 14, 15, 16, 18, 22, 23, 24, 25, 26, 27, 28, 30, 31, 32, 33, 34, 35, 38, 39, 40, 47, 48, 50, 54, 55, 56, 57, 58, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 78, 79, 80, 81, 82, 83, 86, 87, 88, 113, 114, 115, 116, 118, 119, 120, 121, 122, 123, 124, 126, 127, 128, 162, 166, 167, 168, 172, 175, 176, 177, 189, 191, 192, 194, 195, 196, 199, 200, 209, 210, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 260, 261, 262, 263, 264, 266, 267, 268];
%predictors = predictors(:,preSelectedFeatureIdxs);
correlations = corr(predictors(:,preSelectedFeatureIdxs));

exclude = zeros(1,length(preSelectedFeatureIdxs));
numExcludedIndices = 0;

threshold = 0.92;

%exclude indices j for which corr(i,j) > treshold
for i = 1 : length(correlations)
    for j = i + 1 : length(correlations)
        correlation = correlations(i,j);
        if correlation > threshold
            numExcludedIndices = numExcludedIndices + 1;
            excludeIdx = preSelectedFeatureIdxs(j);
            exclude(numExcludedIndices) = excludeIdx;
        end
    end
end

exclude = exclude(1:numExcludedIndices);
exclude = unique(exclude);

preSelectedFeatures = setdiff(preSelectedFeatureIdxs,exclude);

%double check that there is no correlation > threshold
correlations = corr(predictors(:,includeIdxs));
nonDiagIndices = ~eye(size(correlations));
correlations = correlations(nonDiagIndices);
correlations = reshape(correlations,length(nonDiagIndices)-1,length(nonDiagIndices));
maxCorrs = max(correlations,[],1);
plot(1:length(maxCorrs),maxCorrs,'*');
title('correlations');

fprintf('Features after eliminating highly correlated: \n');
for i = 1 : length(preSelectedFeatures)
    feature = preSelectedFeatures(i);
    fprintf('%d, ',feature);
end
%}


%% Wrapper Approach
%run feature selection
preselectedFeatureIdxs = [1, 2, 3, 4, 5, 6, 9, 10, 11, 14, 15, 16, 18, 22, 23, 24, 25, 26, 27, 28, 30, 32, 47, 57, 58, 65, 66, 67, 68, 69, 70, 71, 72, 81, 82, 83, 86, 87, 88, 172, 175, 176, 177, 189, 191, 192, 241, 242, 243, 244, 248, 249, 250, 251, 254, 255, 256, 260, 262, 263, 266, 267];
cv = cvpartition(responses,'kfold',10);

options = statset('display','iter');

%classf = @(xtrain,ytrain,xtest,ytest) sum(ytest ~= classify(xtest,xtrain,ytrain,'quadratic'));
fun = @(XT,yT,Xt,yt) (sum(~strcmp(yt,classify(Xt,XT,yT,'quadratic'))));

fs = sequentialfs(fun,predictors(:,preselectedFeatureIdxs),responses,'cv',cv,'options',options);
preselectedFeatureIdxs(fs)
