# COPI
Classifier of Phosphosite Interactions for predicting phosphorylation sites that stimulate protein-protein interactions

# p-SER and p-THR summary python script
For each entry in the PDB that contains more than 1 distinct protein and at least one p-SER (SEP) or p-THR (TPO). From here, we collected sites at the interface of a PPI, specifically those interacting with exactly one other protein to stay consistent with the binary approach taken by the Predictomes. The two proteins involved in the phospho-interaction must both be at least 10 aa in length to be used with the HMS folding portal. Furthermore, we removed phospho-interactions involving proteins containing non-biological amino acids, as AlphaFold-Multimer cannot predict these. 

# Dedupe pairs python script
Creates two 20 aa windows with 10 aa on either side of the phosphosite and closest binding residue. If the phosphosite and binding windows both have a sequence identity of 80% or greater, they are considered a duplicate and removed from our training set. We used this script to finalize the sites that we eventually folded with AlphaFold-Multimer. 

# Check recovered python script
Checks whether the phosphosite is predicted by AlphaFold-Multimer to lie within 6 Angstroms of the nearest binding residue from the expreimentally determined structure. We used this script to determine "successful" AlphaFold-Multimer predictions of regulatory phosphosites, and these predictions were used as COPI's positive training set

# 
