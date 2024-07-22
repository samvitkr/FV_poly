% © 2022 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
classdef FrequencyVector
    properties (Access = protected)
        frequencyVector
    end
    
    methods
        function obj = FrequencyVector(nSamples, dtSampling)
            arguments
                nSamples (1,1) {mustBeInteger, mustBeNonnegative} = 0
                dtSampling (1,1) double {mustBeReal, mustBeNonnegative} = 0.0
            end
            if (nSamples == 0) || (dtSampling == 0.0)
                % the default values correspond to invalid inputs: the frequency vector is not defined if:
                % (i) there are more than zero samples (nSamples > 0), but the sampling rate is zero (dtSampling == 0), or
                % (ii) there are zero samples (nSamples == 0). In this case the sampling rate is irrelevant
                % in both cases, we set the frequency vector to NaN to indicate the invalid input
                warning('Frequency vector is undefined for given input. Returning NaN.')
                obj.frequencyVector = NaN;
            else
                samplingFrequency = 1 / dtSampling;
                fundamentalFrequency = samplingFrequency / nSamples;
                obj.frequencyVector = (0 : (nSamples-1)) * fundamentalFrequency;  % positive frequencies only
                % convert large positive to negative frequencies
                if mod(nSamples, 2) == 0
                    idxNyquist = nSamples/2 + 1;
                    obj.frequencyVector(idxNyquist : end) = obj.frequencyVector(idxNyquist : end) - samplingFrequency;
                else
                    idxNegStart = (nSamples+1)/2 + 1;
                    obj.frequencyVector(idxNegStart : end) = obj.frequencyVector(idxNegStart : end) - samplingFrequency;
                end
            end
        end
        
        function frequencyVector = get_frequency_vector(obj)
            frequencyVector = obj.frequencyVector;
        end
        
        function radialFrequencyVector = get_radial_frequency_vector(obj)
            radialFrequencyVector = 2.0 * pi * obj.frequencyVector;
        end
        
        function wavelengthVector = get_wavelength_vector(obj)
            radialFrequencyVector = obj.get_radial_frequency_vector();
            wavelengthVector = NaN(size(radialFrequencyVector));
            wavelengthVector(2:end) = (2.0 * pi) ./ radialFrequencyVector(2:end);
        end
    end
end