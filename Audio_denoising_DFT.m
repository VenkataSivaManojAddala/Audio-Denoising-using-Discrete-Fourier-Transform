% Read in the audio file
[s, Fs] = audioread('64_Systolic Split S1 #3.wav');
SNR_input = (10*log10(var(s)/var(randn(size(s)))));
s1 = awgn(s, 10, 'measured'); 
SNR_noisy = (10*log10(var(s1-s)/var(s)));
disp(['Input SNR: ' num2str(SNR_input) ' dB']);
disp(['Noisy SNR: ' num2str(SNR_noisy) ' dB']);

% Define the  required parameters
N = length(s);  % Number of samples in the signal
L = 512;  % Length of the DFT window
K = N/L;  % Number of DFT windows
alpha = 1;  % Weighting factor

% Initialize the output signal
y = zeros(N, 1);

% Perform DFT-based denoising for each window taken
for k = 0:K-1
    % Extracting the current window from the above windows
    x = s1(k*L+1:(k+1)*L);
    
    % Computing the DFT of the above selected window
    X = zeros(L, 1);
    for n = 0:L-1
        for m = 0:L-1
            X(n+1) = X(n+1) + x(m+1)*exp(-1j*2*pi*n*m/L);
        end
    end
    
    % Computeing the magnitude spectrum and then applying the denoising filter
    S = abs(X);
    G = 1./(1 + alpha./S.^2);
    
    % Applying the filter to the DFT of the window
    Y = G.*X;
    
    % Computing the inverse DFT of the filtered signal
    y_window = zeros(L, 1);
    for n = 0:L-1
        for m = 0:L-1
            y_window(n+1) = y_window(n+1) + Y(m+1)*exp(1j*2*pi*n*m/L);
        end
    end
    
    % Add the filtered window to the output signal
    y(k*L+1:(k+1)*L) = y(k*L+1:(k+1)*L) + real(y_window);
end
% Normalizing the output signal
y = y/max(abs(y));

SNR_output = (10*log10(var(y-s)/var(s)));
disp(['Output SNR: ' num2str(SNR_output) ' dB']);

audiowrite('denoised_audio_file.wav', y, Fs);
 subplot(3,1,1);plot(s);ylabel('Amplitude');title('Clean signal')
 subplot(3,1,2);plot(s1);axis tight;ylabel('Amplitude'),title('Noisy signal')
 subplot(3,1,3);plot(y);axis tight;ylabel('Amplitude'),title('Denoised')

