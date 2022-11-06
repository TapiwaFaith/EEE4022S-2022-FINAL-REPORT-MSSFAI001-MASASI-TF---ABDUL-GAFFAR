# EEE4022S-2022-FINAL-REPORT-MSSFAI001-MASASI-TF---ABDUL-GAFFAR

The repository contains the code used to obtain the results in the report.

Step 1: Download the zipped folder that contains the three scripts LoadRadarData.m, computeEntropy.m and computeContrast.m
        LoadRadarData.m is the main program and needs computeEntropy.m and computeContrast.m to be in the same folder in order to run.
        computeEntropy.m is a function used in computing the entropy values.
        computeContrast.m is a function used in computing contrast values.
        
Step 2: Find the following line of code in LoadRadarData.m (line 37)
        HRRProfilesToPlotSubset = HRRProfilesNotAligned(513:640,:);
        Change the range profile values as desired (selecting range profiles to be worked on by the range alignment algorithms).
        
Step 3: Find the following line of code in LoadRadarData.m (line 75)
        r = xcorr(abs(HRRReferenceRangeProfile(129, :)), abs(HRRProfilesToPlotSubset(i, :)));
        Change the range profile number as desired. In this case it has been set to 129 meaning the reference range profile is 129.
        
Step 4: Run the code.
        Once the code has completed running, graphs are going to be displayed.
        The graphs have titles to easily identify which range alignment algorithm resulted in which graph.
        Additionally, the following values are going to be displayed in the Command Window in linear units:
            Entropy_Before
            Contrast_Before
            
            EntropyAfterCrossCorrelation
            ContrastAfterCrossCorrelation
            
            EntropyAfterPeakAlignment
            ContrastAfterPeakAlignment
            
            EntropyAfterHaywood
            ContrastAfterHaywood

Note: The CSIR data used for this work is not included in this public repository.
