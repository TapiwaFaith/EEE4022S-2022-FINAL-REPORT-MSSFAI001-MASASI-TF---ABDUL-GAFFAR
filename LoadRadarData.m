clear all;
close all;

%--------------------------------------------------------------------------------------------------------------------
%DAP_2010-10-12_11-09-28_004_Tugboat_P874_G1_sb_HRR.mat
%DAP_2010-10-14_15-26-43_042_RIB_INBOUND_SPEED_HIGH_P1020_G1_sb_HRR.mat
%DAP_2010-10-14_09-43-33_010_zayaan_inbound_singlebin_P455_G1_sb_HRR.mat

%Loading Data
load('DAP_2010-10-12_11-09-28_004_Tugboat_P874_G1_sb_HRR.mat');             %copy and paste the file containing the data from above
HRRProfilesNotAligned = sb_HRR.G1.HRR_NoMC_calib.';                         %this variable holds all the range profiles
Effective_PRF =  1/(InfoData.TgtData.timePatterns(3)-InfoData.TgtData.timePatterns(2));
RangeResolution_m = sb_HRR.G1.xaxis_downrange_m(3) - sb_HRR.G1.xaxis_downrange_m(2); 

%--------------------------------------------------------------------------------------------------------------------
%Plot the HRR profiles after calibration, but without any motion compensation

y_axis = 1:size(sb_HRR.G1.HRR_NoMC_calib,2);                                %assigns the range profiles to the y-axis for plotting graphs
x_axis = 1:1: size(sb_HRR.G1.HRR_NoMC_calib,1);                             %assigns the range bins to the x-axis for plotting graphs
HRRProfilesToPlot = 20*log10(abs(HRRProfilesNotAligned));                   %scales the unaligned range profiles to decibels values for easy plotting
figure; imagesc(x_axis, y_axis, HRRProfilesToPlot);                         %creates a figure and plots the high range resolution profiles 
set(gca,'FontSize',11)                                                      %sets the font size for the figure
xlabel('Range (bins)','fontsize',11);                                       %labels the x-axis
ylabel('Profile Number','fontsize',11);                                     %labels the y-axis
title('Original HRR profiles','fontsize',11);                               %adds the title to the figure
colormap('jet');                                                            %returns the jet colormap as a three-column array to be put on the colour bar
colorbar;

%--------------------------------------------------------------------------------------------------------------------
%Subset of data to be worked on
HRRProfilesToPlotSubset = HRRProfilesNotAligned(513:640,:);                 %select the subset of range profiles to be worked
HRRProfilesToPlotting = 20*log10(abs(HRRProfilesToPlotSubset));             %scales the unaligned range profiles to decibels values for easy plotting

%--------------------------------------------------------------------------------------------------------------------
%Reference range profile coming from here
HRRReferenceRangeProfile = HRRProfilesNotAligned(1:size(HRRProfilesNotAligned,1),:);

%--------------------------------------------------------------------------------------------------------------------
%Plotting subset of data to be worked on
figure;
imagesc(HRRProfilesToPlotting);                                             %plots the subsets of range profiles to be worked on
hold on;

for i = 1: size(HRRProfilesToPlotSubset, 1)
    [value, index] = max(abs(HRRProfilesToPlotSubset(i, :)));               %adds a black x at the position of the peak scatterer to help in identifying it
    plot(index, i, 'kx', 'MarkerSize',8, 'LineWidth',0.5);
end

set(gca,'FontSize',11)
xlabel('Range (bins)','fontsize',11);
ylabel('Profile Number','fontsize',11);
title('Subset HRR Profiles','fontsize',11);
colormap('jet');
colorbar;

%--------------------------------------------------------------------------------------------------------------------
%Entropy and Contrast Before

%Entropy before range alignment
Entropy_Before = computeEntropy(HRRProfilesToPlotSubset);
Entropy_Before

%Contrast before range alignment
Contrast_Before = computeContrast(HRRProfilesToPlotSubset);
Contrast_Before

%--------------------------------------------------------------------------------------------------------------------
%Cross Correlation

% Align HRR profiles
NumberOfRangeBins = size(HRRProfilesToPlotSubset, 2);

for i = 1: size(HRRProfilesToPlotSubset, 1)
    %find the cross correlation between the reference profile and the profiles to be aligned
    r = xcorr(abs(HRRReferenceRangeProfile(129, :)), abs(HRRProfilesToPlotSubset(i, :)));

    %find the index position of the maximum cross correlation value
    [max_val max_idx]  = max(r);

    %using the index position of the max cross correlation value, finds the number of bin shifts required to align the unaligned profile with the reference profile
    NumberOfBinShiftsRequired = max_idx - NumberOfRangeBins;

    %store the number of bin shifts required in a vector for use in the Haywood algorithm
    NumberOfBinShiftsRequired_Vector(i) = NumberOfBinShiftsRequired;

    %circularly shift the profile to be aligned with respect to the reference profile
    CrossCorrelation(i, :) = circshift(HRRProfilesToPlotSubset(i, : ), NumberOfBinShiftsRequired);
end

% Plot aligned HRR profiles
figure;
imagesc(20*log10(abs(CrossCorrelation)));
colorbar;
xlabel('Range (bins)');
ylabel('Range Profiles');
title('Cross Correlation Aligned Range Profiles');
axis xy;
colormap('jet');

%--------------------------------------------------------------------------------------------------------------------
%Entropy and Contrast After Cross Correlation

%Entropy after range alignment
EntropyAfterCrossCorrelation = computeEntropy(CrossCorrelation);
EntropyAfterCrossCorrelation

%Contrast after range alignment
ContrastAfterCrossCorrelation = computeContrast(CrossCorrelation);
ContrastAfterCrossCorrelation

%--------------------------------------------------------------------------------------------------------------------
%Align HRR profiles using peak Alignment

% 
%for loop to allow iteration through each range profile
for i = 1: size(HRRProfilesToPlotSubset, 1)

    %obtain the peak value and its index number of each range profile
    [value, index] = max(abs(HRRProfilesToPlotSubset(i, :)));

    %find the number of bin shifts needed to shift the peak value to the last index position
    BinShiftsNeeded = -(index - NumberOfRangeBins);

    %circulaly shift the peak value to the last index position
    PeakAlignment(i, :) = circshift(HRRProfilesToPlotSubset(i, :), BinShiftsNeeded);
end

% Plot aligned HRR profiles
figure;
imagesc(20*log10(abs(PeakAlignment)));
hold on;

for i = 1: size(PeakAlignment, 1)
    [value, index] = max(abs(PeakAlignment(i, :)));
    plot(index, i, 'kx', 'MarkerSize',8, 'LineWidth',0.5);
end

colorbar;
xlabel('Range (bins)');
ylabel('Range Profiles');
title('Peak Alignment Aligned Range Profiles');
axis xy;
colormap('jet');

%--------------------------------------------------------------------------------------------------------------------
%Entropy and Contrast After Peak Alignment

%Entropy after range alignment
EntropyAfterPeakAlignment = computeEntropy(PeakAlignment);
EntropyAfterPeakAlignment

%Contrast after range alignment
ContrastAfterPeakAlignment = computeContrast(PeakAlignment);
ContrastAfterPeakAlignment

%--------------------------------------------------------------------------------------------------------------------
%Sub-Integer Haywood

%number of elements. NumberOfBinShiftsRequiredVector was obtained from the cross correlation method.
NumberOfElements = length(NumberOfBinShiftsRequired_Vector);

%profile vector contains the number of elements
ProfileVector = 1:1:NumberOfElements;

%low order polynomial value
Order = 1;

%fitting a low order polynomial on the data
Coefficients = polyfit(ProfileVector, NumberOfBinShiftsRequired_Vector, Order);

%finding the non integer bin shifts
b = polyval(Coefficients, ProfileVector);

n = size(HRRProfilesToPlotSubset,2);
m = 0:1:(n-1);

for pos = 1: size(b, 2)
    s = b(1,pos);
    
    %PhaseShiftVector contains the matrix of sub-integer range bin values to shift the signal
    PhaseShiftVector = exp(-1i*2*pi*b(1,pos)*m/n);

    %Applying Fast Fourier Transform to unaligned profiles
    Multiplier = fft(HRRProfilesToPlotSubset(pos, :));

    %Multiplication by a phase ramp which is the phase shift vector
    HaywoodShiftedProfile(pos, :) = ifft(Multiplier.*PhaseShiftVector);
end

% Plot aligned HRR profiles
figure;
imagesc(20*log10(abs(HaywoodShiftedProfile)));
colorbar;
xlabel('Range (bins)');
ylabel('Range Profiles');
title('Haywood Aligned Range Profiles');
axis xy;
colormap('jet');

%--------------------------------------------------------------------------------------------------------------------
%Entropy and Contrast After Haywood

%Entropy after range alignment
EntropyAfterHaywood = computeEntropy(HaywoodShiftedProfile);
EntropyAfterHaywood

%Contrast after range alignment
ContrastAfterHaywood = computeContrast(HaywoodShiftedProfile);
ContrastAfterHaywood

%--------------------------------------------------------------------------------------------------------------------
%Writing elements to the table in excel

T = table(Entropy_Before, EntropyAfterCrossCorrelation, EntropyAfterPeakAlignment, EntropyAfterHaywood, Contrast_Before, ContrastAfterCrossCorrelation, ContrastAfterPeakAlignment, ContrastAfterHaywood);
%writetable(T, 'Tugboat Experimental Results.xlsx', 'WriteMode', 'append', 'WriteVariableNames', false, 'WriteRowNames', true);

%--------------------------------------------------------------------------------------------------------------------
%Compares the different range profiles  to check whether they have change or not the amplitudes of the range profiles

ReferenceRangeProfile = abs(HRRReferenceRangeProfile(10,:));
RangeProfile = abs(HRRReferenceRangeProfile(120,:));
RangeProfile2 = abs(HRRReferenceRangeProfile(769,:));
  
figure;
set(gca,'FontSize',20);
plot(x_axis, ReferenceRangeProfile, 'r-');
hold on;
plot(x_axis, RangeProfile, 'b-');
hold on;
plot(x_axis, RangeProfile2, 'c-');
xlabel('Range (bins)','fontsize',11);
ylabel('Amplitude','fontsize',11);
legend('Range profile 10', 'Range profile 120', 'Range profile 769', 'fontsize',11);