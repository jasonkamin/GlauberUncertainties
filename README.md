# GlauberUncertainties
Set of macros to run over TGlauberMC simulations (ROOT TTree format) and calculate relative uncertainties for NColl, NPart, and TAA. Typically, 1M events is the right ballpark for reasonable statistics.  For the density parameter variations, we've been using 100 variations.  This number could be bumped up but then everything would take more time and 100 appears to be enough for robust results. 

We calculate the relative uncertainties for NColl, NPart, and TAA. 
We can obviously calculate the central values as well but, currently, 
these values are all based on an impact-parameter slicing of the centralities
whereas, for CMS, we actually use an HF-response based slicing. Therefore, 
we use the relative uncertainties and use the central values from the HF-based method. 

The current procedure is to calculate the uncertainties based on 4 "sources": 

1. Nucleon Density Parameters: 
   varying the p&n density parameters (by gaussians with 1 sigma; 100 variations) and looking at how the variables (statistically) respond. 
2. Minimum Imposed Distance between Nucleons. 
   varied by a factor of 2 up and down. 
3. Inelastic Scattering Cross-Section
   vary the sigma\_inel by ~5%, up and down. 
4. HF efficiency and Electromagnetic Contamination. 
   vary efficiency + EM contamination by 2%, up and down. 

The procedure is to run MakeLatexTables.C with the appropriate list of vectors (which define what variations we'll look at) for 2,3,4. 
For 1, run Plot_pnPars_uncert.C with the appropriate flags set. 
