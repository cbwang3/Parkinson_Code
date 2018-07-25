classdef FeatureExtractor < handle
    
    %TODO  entropy, spectralFlatness
    
    properties (Access = public)
        nFeatures;
        featureNames;
    end
    
    properties (Access = private)
        featureExtractors = {@min,@max,@mean,@var,@std,@median,@trapz,@aav,...
            @mad,@iqr,@skewness,@kurtosis,@rms};%TODO add energy here
        
        %fourierFeatureExtractors = {@maximumFrequency, @spectralEnergy, @spectralEntropy,@spectralCentroid,@spectralSpread,...
            %@dc};
    end
    
    methods (Access = public)
        
        function obj = FeatureExtractor()
           segment = rand(50,7);
            [~,obj.featureNames] = obj.extractFeaturesForSegment(segment);
            obj.nFeatures = length(obj.featureNames);
        end
        
        function [featureVector, featureNames] = extractFeaturesForSegment(obj,segment)
            %% Calculate additional signals.
            %signal columns are: 1: timestamp , 2: accelerationX, 3: accelerationY, 4:
            %accelerationZ, 5: rotationX, 6: rotationY, 7: rotationZ, 8:
            %accelerationMagnitude, 9: rotationMagnitude
            
                signalNames = {'Ax','Ay','Az','Gx','Gy','Gz','MA','MG'};
            
            accelerationMagnitude = single(segment(:,2).^2 + segment(:,3).^2 + segment(:,4).^2);
            rotationMagnitude = single(segment(:,5).^2 + segment(:,6).^2 + segment(:,7).^2);
            segment = [segment, accelerationMagnitude];
            segment = [segment, rotationMagnitude];
            clear accelerationMagnitude;
            clear rotationMagnitude;
            
            %fourierTransform = dsp.FFT('FFTImplementation','Radix-2','Normalize',false,'FFTLengthSource','Property','FFTLength',256,'WrapInput',false);
            
            %% Extract features               
            featureVector = zeros(1,obj.nFeatures);
            featureNames = cell(1,obj.nFeatures);
            
            featureCounter = 1;
            for currentFeature = 1 : length(obj.featureExtractors)
                for currentSignal = 2 : size(segment,2)
                    featureExtractorHandle = obj.featureExtractors(currentFeature);
                    featureExtractorHandle = featureExtractorHandle{1};
                    featureVector(featureCounter) = featureExtractorHandle(double(segment(:,currentSignal)));
                    featureNames(featureCounter) = obj.convertFeatureToString(featureExtractorHandle,currentSignal,signalNames);
                    featureCounter = featureCounter + 1;
                end
            end
            
            %quantile 4 * nSignal features = 24
            numQuantileParts = 4;
            for currentSignal = 2 : size(segment,2)
                quantilesResult = quantile(segment(:,currentSignal),numQuantileParts);
                featureStringName = obj.convertFeatureToString(@quantile,currentSignal, signalNames);
                
                for quantilePart = 1 : numQuantileParts
                    featureVector(featureCounter) = quantilesResult(quantilePart);
                    finalFeatureName = sprintf('%s%d',featureStringName{1},quantilePart);
                    featureNames(featureCounter) = {finalFeatureName};
                    featureCounter = featureCounter + 1;
                end
            end
            
            % zero crossing: only for ax,ay,az
            for currentSignal = 2 : 4
                featureVector(featureCounter) = zrc(segment(:,currentSignal));
                featureNames(featureCounter) = obj.convertFeatureToString(@zrc,currentSignal,signalNames);
                featureCounter = featureCounter + 1;
            end
            
            % octants: 2
            accelerationOctant = octant(segment(:,2), segment(:,3), segment(:,4));
            rotationOctant = octant(segment(:,5), segment(:,6), segment(:,7));
            
            [~, maxAccelIndex] = max(segment(:,8));
            featureVector(featureCounter) = accelerationOctant(maxAccelIndex);
            featureNames(featureCounter) = {'accelMagnitudeOctant'};
            featureCounter = featureCounter + 1;
            [~, maxRotationIndex] = max(segment(:,9));
            featureVector(featureCounter) = rotationOctant(maxRotationIndex);
            featureNames(featureCounter) = {'rotationMagnitudeOctant'};
            featureCounter = featureCounter + 1;
            
            %sma acceleration: 1
            featureVector(featureCounter) = sma(segment(:,2:4));
            featureNames(featureCounter) = {'smaAcceleration'};
            featureCounter = featureCounter + 1;
            
            %sma rotation: 1
            featureVector(featureCounter) = sma(segment(:,5:7));
            featureNames(featureCounter) = {'smaRotation'};
            featureCounter = featureCounter + 1;
            
            %svm acceleration: 1
            featureVector(featureCounter) = svmFeature(segment(:,8));
            featureNames(featureCounter) = {'svmAcceleration'};
            featureCounter = featureCounter + 1;
            
            %svm rotation: 1
            featureVector(featureCounter) = svmFeature(segment(:,9));
            featureNames(featureCounter) = {'svmRotation'};
            featureCounter = featureCounter + 1;
            
            %correlation coefficients: 8 * 7 = 56
            corrCoeffResult = corrcoef(segment(:,2:end));
            for i = 1 : size(corrCoeffResult,1)
                for j = i + 1 : size(corrCoeffResult,2)
                    featureVector(featureCounter) = corrCoeffResult(i,j);
                    featureNames(featureCounter) = obj.getCorrelationFeatureString(i,j,'corr',signalNames);
                    featureCounter = featureCounter + 1;
                end
            end
            
            %cross correlation coefficients: 8 * 7
            for i = 2 : size(segment,2)
                for j = i + 1 : size(segment,2)%for each pair of signals
                    featureVector(featureCounter) = maxCrossCorr(segment(:,i),segment(:,j));
                    featureNames(featureCounter) = obj.getCorrelationFeatureString(i-1,j-1,'xcorr',signalNames);
                    featureCounter = featureCounter + 1;
                end
            end
            
            %{
            %fourier features
            for currentFeature = 1 : length(obj.fourierFeatureExtractors)
                for currentSignal = 2 : size(segment,2)
                    featureExtractorHandle = obj.fourierFeatureExtractors{currentFeature};
                    window = complex(segment(:,currentSignal));
                    featureVector(featureCounter) = featureExtractorHandle(window,fourierTransform);
                    featureNames(featureCounter) = obj.convertFeatureToString(featureExtractorHandle,currentSignal,signalNames);
                    featureCounter = featureCounter + 1;
                end
            end
            %}
            %disp(featureCounter);
            % duration
            %featureVector(featureCounter) = size(segment,1);
            %featureNames(featureCounter) = {'duration'};
        end
    end
    
    methods (Access = private)
        function [featureString] = convertFeatureToString(~,featureExtractorHandle,signalIdx, signalNames)
            signalName = signalNames(signalIdx-1);
            featureString = sprintf('%s%s',func2str(featureExtractorHandle),signalName{1});
            featureString = {featureString};
        end
        
        function [featureString] = getCorrelationFeatureString(~,row,col, correlationType, signalNames)
            signalName1 = signalNames(row);
            signalName2 = signalNames(col);
            featureString = sprintf('%s%s%s',correlationType,signalName1{1},signalName2{1});
            featureString = {featureString};
        end
    end
end


