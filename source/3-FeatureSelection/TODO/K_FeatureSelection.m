close all;

path = genpath('./');
addpath(path);
loadClasses;

fileNames = {'1-lukas','2-dino','3-luis', '4-mauricio','5-leo','7-noa','8-alex','9-jakob','10-michi','11-roko','12-theo','13-aldin'};
%fileNames = {'2-dino'};

%feature selection for segmentation
global exercisesTable;
lazyLoadExercisesTableForSegmentation(fileNames);

%feature selection for classification
%labelingStrategy = GroupSimilarClassesLabelingStrategy();
%labelingStrategy = SprintJoggingExerciseLabelingStrategy();
%segmentationAlgorithm = ManualSegmentation();   
%exercisesTable = createTableForClassification(fileNames,segmentationAlgorithm,labelingStrategy);

%prepare data
nFeatures = size(exercisesTable,2);
predictors = exercisesTable(:,1:nFeatures-1);
predictors = table2array(predictors);
responses = exercisesTable.label;

%run feature selection
mdl = fscnca(predictors,responses);

%plot results
figure();
plot(mdl.FeatureWeights,'ro');
grid on
xlabel('Feature index');
ylabel('Feature weight');

%print results
[sortedWeights, sortedFeatures] = sort(mdl.FeatureWeights,'descend');
nFeatures = 50;

fprintf('Best features are:\n');

for i = 1 : nFeatures
    featureIdx = sortedFeatures(i);
    featureName = exercisesTable.Properties.VariableNames(featureIdx);
    weight = sortedWeights(i);
    fprintf('%d- (%d-%s) - %.2f\n',i, featureIdx,featureName{1},weight*10000);
end

%print feature list to copy into your code
for i = 1 : nFeatures
    featureIdx = sortedFeatures(i);
    fprintf('%d, ',featureIdx);
end
