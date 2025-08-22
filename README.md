# Temperature-based-spectral-indices
We present our MatLab code for calculating temperature-eigenvalues-based topological indices.
Here is the flowchart:

Step 1: Use HyperChem drawing module to construct 3D molecular graph of $\Phi$. It delivers a file with .hin extension.

Step 2: Feed .hin file to TopoCluj to compute the distance matrix of $\Phi$ and generate .m file corresponding to .hin file.

Step 3: Compute all temperature-related spectral descriptors by inputting .m to MatLab and employ our code written in MatLab. 
