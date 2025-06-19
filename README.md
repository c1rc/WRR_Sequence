# WRR_Sequence
CVFEM groundwater flow model modified to work with Sequence

NSF_CVFEM_Sequence.m is a modified version of the CVFEM dynamic grid groundwater flow model published in Voller (2025). This version is tailored to run with Sequence (https://github.com/sequence-dev/sequence) outputs. This is the only file that requires input from the user to run this linked version of the model. The other two files are for two separate purposes-
a)sequence_import_3_generalized.m - Imports top and bottom boundary from Sequence
b)plotstfac_to_voller.m - Imports the stratigraphy from facies as polygon shapes into the main code
