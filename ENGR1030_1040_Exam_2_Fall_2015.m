%% Fall 2015 ENGR 1030/1040 Exam 2  (100 Points)
% Michael Davenport
% Student's Name: ________________________  
% Student's Section:
% Date: 

% General Instructions
% Many years ago I took a course on wavelets, which are simple oscillating
% functions of a finite duration. The wavelet transform applies wavelets to
% a wide variety of applications in digital signal processing. For my
% course project I experimented with using wavelets as a frequency filter
% on an audio signal like bands of an equalizer. In the original MATLAB
% file that I wrote for the project I changed the various parameters used
% in the program manually and re-ran it each time I wanted to experiment
% with different settings.
%
% For this test you will use the user-controlled input/output functions
% from Chapter 7, and the selection and repetition structures from Chapters
% 8 and 9 to provide an interactive command line interface for automating
% the testing process from the MATLAB command prompt. The instructions for
% each part of the test are given at the top of the cell. I've also added
% comments within the cell where your code modifications should be placed.

%% Part 1: Selecting an audio file and downsampling
% a) (5 Points) Add code to allow the user to enter the name of an audio
% file so that the filename is not hardcoded in the audioread() function.
% For this test you can use the attached file 'CanYouFeelIt.wav')

% b) (10 Points) Add code to allow the user to enter an optional downsampling
% value. The program should downsample the original audio file if the user
% enters an intger greater than 1. (Hint: round up the user input so that
% it will be an integer.) Since you may be rounding up the number, be sure
% to output an message telling the user what the actual downsampling value
% will be, but only if downsampling actually occurs. Also, be sure to
% update the title of the audio sequence plot to indicate if the audio data
% is the original or downsampled.
%
% c) (20 Points) Add a test on the user input value for downsampling to make
% sure it is a number. (Hint: use the MATLAB built-in function isnumeric to
% perform the test.) If the user does not enter a number, output an error
% message and ask them to re-enter the number. (Hint: you will need to use
% a repetition structure that will continue to ask for user input until a
% valid number is received.)
clear, clc

%%% Add your user input code for Part 1 in this section

%%% end
runfile = 5;
%%% Make any necessary modifications to the code in this section for Part 1
while runfile ~=1;
% Load the audio file and extract data
fileTypes = {'*.mp3';'*.wav';'*.wma';'*.ogg';'*.flac';'*.au';'*.aiff';'*.aif';
    '*.aifc';'*.mp4';'*.m4a';};
filename = uigetfile(fileTypes,'Pick an AUDIO file');



[audio_signal, sampling_frequency] = audioread(filename);


%play file original file if wanted 
%gong = audioplayer(audio_signal, sampling_frequency);
%play(gong);
[length_audio,channels]=size(audio_signal);

leftChannel = audio_signal(:,1);
if channels > 1;
rightChannel = audio_signal(:,2);
else rightChannel = audio_signal(:,1);
end
x_len = length(audio_signal);



% Downsample the original input
b = menu('would you like to downsample','yes','no');

if b == 1
    c = menu('What factor would you like to downsample by','1','2','3');
    downsampling_factor = c;
	%downsampling_factor = 2;
    leftChannel = leftChannel(1:downsampling_factor:x_len);
    rightChannel = rightChannel(1:downsampling_factor:x_len);
    audio_signal = [leftChannel'; rightChannel'];
    audio_signal = audio_signal';
    x_len = (x_len - 1)/ downsampling_factor;					
	sampling_frequency = sampling_frequency / downsampling_factor;												%sampling frequency of downsampled sequence
end

% plot the audio sequence
figure
subplot(2,1,1)
plot (leftChannel)
title('Original Left Channel')
xlabel('Number of Samples')
ylabel('Signal Magnitude')
text(0,0, ['Sampling frequency: ' num2str(sampling_frequency) ' Hz'])
subplot(2,1,2)
plot (rightChannel)
title('Original Right Channel')
xlabel('Number of Samples')
ylabel('Signal Magnitude')
text(0,0, ['Sampling frequency: ' num2str(sampling_frequency) ' Hz'])

%%% End of code affected for Part 1

%% Part 2: Optional FFT Plot
% a) (5 Points) Give the user the option of displaying the FFT plot of the
% original signal. If the user enters either an upper case Y or lower case
% y, then plot the FFT data.

% b) (10 Points) If the user has selected to display the FFT plot, ask them
% to retrieve the frequency with the largest magnitude by clicking on the
% highest peak of the graphs of each channel. They should only be allowed
% to click two times, once per graph. Display the results in a statement
% indicating the frequencies with the largest magnitude.

%%% Add your user input code for Part 2 in this section

%%% end

% Compute the Fast Fourier Transform (FFT)
% This section calculates the frequency content of each channel of the
% audio signal.
% Calculate the number of points for the FFT. We want a value that is a
% power of 2 and just larger than the number of samples in the audio
% signal.
N = 2^nextpow2(x_len);

% Calculate the frequency axis of the N-point FFT
faxis = sampling_frequency*(0:(N/2-1))/N;

%Generate FFT of the input signal
rightChannel = fft(rightChannel, N);
leftChannel = fft(leftChannel, N);

%%% Make any necessary modifications to the code in this section for Part 2

% Plot the FFT data
if 1
   r_mag = abs(rightChannel);
   figure
   subplot(2,1,1)
   semilogx(faxis, r_mag(1:(N/2)))
   title('Original Right Channel Frequency Content')
   xlabel('Frequency (Hz)')
   grid
   l_mag = abs(leftChannel);
   subplot(2,1,2)
   semilogx(faxis, l_mag(1:(N/2)))
   title('Original Left Channel Frequency Content')
   xlabel('Frequency (Hz)')
   grid
end

%%% End of code affected for Part 2

%% Part 3: Selecting a Wavelet for the Transformation
% a) (25 Points) Let the user choose a wavelet function to use for the
% transformation. If they enter a '1', use the Haar Wavelet; if '2', use
% the Mexican Hat Wavelet; and if '3', use the Morlet Wavelet. Use the
% switch-case structure to perform the selection.

% b) (5 Points) Print out the name of the wavelet being used for the
% transformation and also the calculated value of the center_frequency
% variable. Also, update the title of the filtered audio plot at the end to
% indicate whether it is using the original file or a downsampled version.

% c) (10 Points) Use the variable where you previously stored the user input
% to determine whether or not to display the FFT plot of the original
% signal to again determine if the FFT plot of the transformed signal
% should be displayed.

% Extra Credit d) (10 Points) Modify the code from a) to allow the user to
% input a value of 0 to repeat the wavelet transform three times using all
% three wavelets. You will need to use a repetition structure around the
% switch-case structure.

%%% Add your user input code for Part 3 in this section


% GUI menu to pick a wavelet to use

Wavepick = menu('Choose a wavelet to apply ','Haar Wavelet',...
'Mexican Hat Wavelet','Morlet Wavelet','Poisson Wavelet 4', 'Real Shannon Wavelet' );

%Rich Gordon -- add an option to show which graphs you'd like to see
show_wavelet = input('would you like to see the original wavelet waveform? y/n');
show_wavelet_dft = input('would you like to see the original wavelet dft? y/n');


switch Wavepick
    
    case 1
        %Haar Wavelet
        option = 1;
        
    case 2
        %Mexican Hat Wavelet
        option = 2;
        
    case 3
        %Morlet Wavelet
        option = 3;
        
    case 4 
        %Wavelet Number 4
        option = 4;
        
   case 5
        % No wavelet
        option = 5;
    case 0
        % dialog box is exited without making a choice
        % default to Morlet Wavelet
        option = 3;
   casozero = msgbox('No choice made. Defaulting to Morlet Wavelet');
    otherwise
        
        disp('You broke my code. How inconsiderate!')
end
    


%%% end

%%% Make any necessary modifications to the code in this section for Part 3

% Compute the wavelet function, determine its FFT and center frequency
a = 1;			% dialation factor, 1 = mother wavelet


switch option
	case 1	%Haar Wavelet
        t = 0:(1/sampling_frequency):1;
		haar = abs(a)/2 - (0:(1/sampling_frequency):abs(a));
		haar = sign(a)*sign(haar)/sqrt(abs(a));
        if show_wavelet == 'y'
            figure
            plot(t,haar)
            grid
            title('Haar Mother Wavelet')
        end
        haar_fft = fft(haar, N);
        haar_mag = abs(haar_fft);
        if show_wavelet_dft == 'y'
            figure
            semilogx(faxis, haar_mag(1:(N/2)))
            title('DFT of Haar Mother Wavelet')
            xlabel('Frequency (Hz)')
            grid
        end
        [val, int] = max(haar_mag(1:(N/2)));
        center_frequency = faxis(int);   
        w_fft = haar_fft';
        w_fft = w_fft / max(w_fft);
	case 2	%Mexican Hat Wavelet
		t = 0:(1/sampling_frequency):10;		%time axis for the wavelet function
        mex_hat = 1/sqrt(a)*(1 - 2*(5*t/a).^2).*exp(-(5*t/a).^2);
        if show_wavelet == 'y'		%1 = show plot of mother wavelet
            figure
            plot(t,mex_hat)
            grid
            title('Mexican Hat Mother Wavelet')
        end
        mex_fft = fft(mex_hat, N);
        mex_mag = abs(mex_fft);
        if show_wavelet_dft == 'y'
            figure
            semilogx(faxis, mex_mag(1:(N/2)))
            title('DFT of Mexican Hat Mother Wavelet')
            xlabel('Frequency (Hz)')
            grid
        end
        [val, int] = max(mex_mag(1:(N/2)));
        center_frequency = faxis(int);
        w_fft = mex_fft';
        w_fft = w_fft / max(w_fft);
   case 3	%Real Valued Morlet Wavelet	
        a = 0.88319778442383/15000;	%used for generating wavelet functions
        t = 0:(1/sampling_frequency):10;
        morlet = 1/sqrt(a)*exp(-0.754*((t/a).*2)).*cos((pi*sqrt(2/log(2))).*(t/a));
        if show_wavelet == 'y'
            figure
            plot(t,morlet)
            grid
            title('Morlet Mother Wavelet')
        end
        mor_fft = fft(morlet, N);
        mor_mag = abs(mor_fft);		%normalize the sequence
        if show_wavelet_dft == 'y'
            figure
            semilogx(faxis, mor_mag(1:(N/2)))
            title('DFT of Morlet Mother Wavelet')
            xlabel('Frequency (Hz)')
            grid
        end
        if show_wavelet_dft == 'y'
            mor_fft1 = mor_fft / max(mor_fft);
            mor_mag = abs(mor_fft1);
            figure
            semilogx(faxis, mor_mag(1:(N/2)))
            title('Normalized DFT of Modified Morlet Wavelet')
            xlabel('Frequency (Hz)')
            grid
        end
        [val, int] = max(mor_mag(1:(N/2)));
        center_frequency = faxis(int);
        w_fft = mor_fft';
        w_fft = w_fft / max(w_fft);
        case 4 %Poisson Wavelet -- Rich Baird
        %M(t)= 1/pi * (1-t^2)/(1+t^2)^2
        t = 0:(1/sampling_frequency):1;
		% apply function to every t
        poisson = 1/pi .* (1-t.^2) ./ (1+t.^2).^2;
        if show_wavelet == 'y'
            figure
            plot(t,poisson)
            grid
            title('Poisson Mother Wavelet')
        end
        poisson_fft = fft(poisson, N);
        poisson_mag = abs(poisson_fft);
        if show_wavelet_dft == 'y'
            figure
            semilogx(faxis, poisson_mag(1:(N/2)))
            title('DFT of Poisson Mother Wavelet')
            xlabel('Frequency (Hz)')
            grid
        end
        [val, int] = max(poisson_mag(1:(N/2)));
        center_frequency = faxis(int);   
        w_fft = poisson_fft';
        w_fft = w_fft / max(w_fft);
    case 5 %Real Shannon Wavelet -- Rich Baird
        %M(t)= sinc(t/2) * cos(3pit/2)
        t = 0:(1/sampling_frequency):1;
		% apply function to every t
        sha = sinc(t./2) .* cos((3 .* pi .* t)./2); 
        if show_wavelet == 'y'
            figure
            plot(t,sha)
            grid
            title('Real Shannon Mother Wavelet')
        end
        sha_fft = fft(sha, N);
        sha_mag = abs(sha_fft);
        if show_wavelet_dft == 'y'
            figure
            semilogx(faxis, sha_mag(1:(N/2)))
            title('DFT of Real Shannon Mother Wavelet')
            xlabel('Frequency (Hz)')
            grid
        end
        [val, int] = max(sha_mag(1:(N/2)));
        center_frequency = faxis(int);   
        w_fft = sha_fft';
        w_fft = w_fft / max(w_fft);
    otherwise
end

%%% Add code here to display which wavelet was selected and the value of
%%% the center_frequency variable


%%% end

% Perform the wavelet transformation
filteredLeftChannel = leftChannel .* w_fft;
filteredRightChannel = rightChannel .* w_fft;

if 1   
    l_mag = abs(filteredLeftChannel);
    figure
    subplot(2,1,1)
    semilogx(faxis, l_mag(1:(N/2)))
    title('Filtered Left Channel Frequency Content')
    xlabel('Frequency (Hz)')
    grid
    r_mag = abs(filteredRightChannel);
    subplot(2,1,2)
    semilogx(faxis, r_mag(1:(N/2)))
    title('Filtered Right Channel Frequency Content')
    xlabel('Frequency (Hz)')
    grid
end

% Finally, convert back to time-domain
newLeftChannel = ifft(filteredLeftChannel, N);
newRightChannel = ifft(filteredRightChannel, N);
filtered_audio = [real(newLeftChannel) real(newRightChannel)];
[length_audio,channels]=size(audio_signal);
if channels == 2;
filtered_audioplay = (max(audio_signal)/max(filtered_audio))*filtered_audio;
else if channels == 1
 filtered_audioplay = (max(audio_signal)/max(filtered_audio(:,1)))*filtered_audio;   
    else
        display('what kind of audio file is this I did not code for this')
    end
end
filtered = audioplayer(filtered_audioplay, sampling_frequency);
play(filtered);

% So we came up with two different ways to convert and play the audio   
% comment 'play(filtered)' above and uncomment the three lines of code below 
% to hear the other method (yes they turned out different) I'm not sure 
% which one is more accurate but I assume it's not mine, it seems too quite
% for the Haar wavelet.


% Write filtered audio signal back to an audio file then play


   % Convert to audio file

%audiowrite('filtered_audio.wav',filtered_audio, sampling_frequency)

   % read converted file

%[soundplay,spr] = audioread('filtered_audio.wav');

   % Play sound 

%sound(soundplay,spr);



% plot the filtered audio sequence for comparison to the original audio
figure
subplot(2,1,1)
newLeftChannel = filtered_audio(:,1);
plot (newLeftChannel)
title('Filtered Original Left Channel')
xlabel('Number of Samples')
ylabel('Signal Magnitude')
subplot(2,1,2)
newRightChannel = filtered_audio(:,2);
plot (newRightChannel)
title('Filtered Original Right Channel')
xlabel('Number of Samples')
ylabel('Signal Magnitude')

 runfile = menu('would you like to run the program again','no','yes');

end
%% Part 4: Matrix Algebra

% a) (10 points) Using the linear equations in the handout containing the
% electrical circuit, use matrix algebra techniques and MATLAB to solve for
% the unknown voltages A and B.

