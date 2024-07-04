#README for LipidIN

###Introduction
The LipidIN is an automated tool designed to facilitate the rapid and accurate identification of lipids in biological samples. This software is applicable to a variety of instruments and uses advanced models, algorithms, and AI technology to automatically analyze the composition of lipids in samples. By using this software, users can identify lipids without the need to establish a standard mass spectrometry library, saving significant time and effort while improving the accuracy and reliability of lipid identification.

###Features
Mass Spectrometry Peak Processing Module: Utilizes existing R packages to process mass spectrometry data in mzML format.
Non-Prior Information Secondary Matching Module: Performs secondary matching with theoretical or real mass spectrometry libraries and normalizes the matching results.
Heuristic Search Module: Based on the relative position of primary information, it conducts heuristic searches using secondary matching scores as prior information to re-evaluate high-score matches.
Prior Information-Based Secondary Matching Module: Enhances poor secondary matching results using primary relative position relationships for better accuracy.
Reverse Lipid Fingerprint Spectrogram Module：Predict  lipid fingerprint spectrogram using the model  we designed inspired by KAN and Muit-head attention.

###System Architecture
The system  main development languages being R，python and C++. While R ,python handles backend processes and C++ accelerates the program. The software employs greedy secondary matching algorithms, heuristic search algorithms, and prior information-based spectrum enhancement algorithms for mass spectrometry analysis, outputting the identification results in CSV format.Using the model to generate reverse lipid fingerprint spectrograms, independent of sample matrices, instruments.

###Modules Description
Mass Spectrometry Peak Processing Module: Processes mzML format data, centralizes peaks, and converts them into list format for easier use.
Non-Prior Information Secondary Matching Module: Matches mass spectrometry peaks with standard libraries using primary and secondary information.
Heuristic Search Module: Reassesses high-confidence matches based on primary information relationships.
Prior Information-Based Secondary Matching Module: Enhances and rematches poor spectrum peaks based on primary relative positions.
Reverse Lipid Fingerprint Spectrogram Module:Generates the reverse lipid fingerprint spectrograms.

###Usage Instructions
Data Input: Select and upload the mzML format file.
Parameter Selection: Choose appropriate parameters for peak processing.
Secondary Matching: Upload the secondary database and set matching parameters.
Heuristic Search: Perform heuristic search based on secondary matching results.
Enhanced Matching: Re-evaluate low-confidence matches with prior information and output final results.
Reverse Lipid Fingerprint Spectrogram Predicting:Using the model to generate reverse lipid fingerprint spectrograms, independent of sample matrices, instruments.
