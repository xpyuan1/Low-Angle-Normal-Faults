# Low-Angle-Normal-Faults
# The codes are used to reproduce Figure 2 and Figure 5 in "Yuan X.P., Olive J.-A., Braun J., Partially-locked low-angle normal faults in cohesive upper crust, Tectonics, 2020."

To construct plots of critical dip as a function of LANF friction (Figure 2), one must run the code with input file BetaMuF.in 
The following command will compile the code: 
gfortran -O BetaMuF.f -o BetaMuF
And this command will run the code with its corresponding input file: 
./BetaMuF<BetaMuF.in 

To construct plots of fault cohesion as a function of fault friction (Figure 5), one must use input file CFMuF.in 
The following commands compile and run the code: 
gfortran -O CFMuF.f -o CFMuF 
./CFMuF<CFMuF.in 

More details on the parameters can be found in the Supplementary Information of Yuan et al. (2020, Tectonics).
