# WRR_Sequence
CVFEM groundwater flow model modified to work with Sequence

NSF_CVFEM_Sequence.m is a modified version of the CVFEM dynamic grid groundwater flow model published in Voller (2025). This version is tailored to run with Sequence (https://github.com/sequence-dev/sequence) outputs. This is the only file that requires input from the user to run this linked version of the model. The other two files are for two separate purposes-
a)sequenceImport.m - Imports top and bottom boundary from Sequence
b)plotStratigraphy.m - Imports the stratigraphy from facies as polygon shapes into the main code

# Julia-Version Branch
NSF_CVFEM_Sequence.jl is the Julia version of the original code with newly added depth varying porosity for aquifers and confining units. This code does not require additional files but it does require the libraries to be installed in Julia.
