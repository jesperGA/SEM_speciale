function [amplitude, frequency] = fourierTransform(force, time)
    % fourierTransform: Compute the amplitude spectrum of a signal using the Fourier transform.
    %
    % Input:
    %   force: Input signal (time-domain).
    %   time: Time vector corresponding to the signal.
    %
    % Output:
    %   amplitude: Amplitude spectrum of the input signal.
    %   frequency: Frequency vector corresponding to the amplitude spectrum.

    % Compute the Fourier transform of the input signal
    fft_result = fft(force) / length(time);
    
    % Calculate the sampling frequency
    sampling_frequency = 1 / (time(2) - time(1));
    
    % Number of samples in the signal
    num_samples = numel(force);
    
    % Generate the frequency vector
    frequency = (0:num_samples - 1) * (sampling_frequency / num_samples);
    
    % Keep only the positive frequencies (up to Nyquist frequency)
    frequency = frequency(1:floor(end / 2));
    
    % Keep only the corresponding part of the amplitude spectrum
    amplitude = abs(fft_result(1:floor(end / 2), :));
end