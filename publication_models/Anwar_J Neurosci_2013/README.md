Stochastic calcium spikes models explanatory note


Model Published in:

Anwar H, Hepburn I, Nedelescu H, Chen W, De Schutter E (2013) Stochastic calcium mechanisms cause dendritic calcium spike variability. J Neurosci (in press)

Scripts Authors: 

Haroon Anwar: haroon.anwar@gmail.com
Iain Hepburn: ihepburn@oist.jp

Software Requirements:

These scripts run on STEPS 2.x, which can be found at http://steps.sourceforge.net.

List of Meshes:

•	5 cylindrical meshes of lengths 10 µm, 20 µm, 40 µm, 80 µm and 160 µm can be found in folder named 'meshes'. These cylindrical meshes can be used with scripts 1-9.
•	A dendritic mesh (‘chop_mesh.inp’) can also be found in folder named "meshes". This mesh can be used with scripts 10 and 11.  

List of Scripts:

1.	The script ‘StochasticCaburst.py’ was used with all cylindrical meshes mentioned above to simulate calcium bursts shown in Figure 2A-C. The traces generated using this script were used for analysis shown in Figure 2 E, Figure 4 and Figure 8D-F.
2.	The script ‘StochasticHH.py’ was used with all cylindrical meshes mentioned above to simulate HH spikes shown in Figure 2D. The traces generated using this script were used for analysis shown in Figure 2F, Figure 3.
3.	The script ‘HybridCaburst_stochCaP.py’ was used cylindrical mesh of 10 µm length to simulate calcium bursts shown in Figure 5A.
4.	The script ‘HybridCaburst_stochCaT.py’ was used cylindrical mesh of 10 µm length to simulate calcium bursts shown in Figure 5B.
5.	The script ‘HybridCaburst_stochBK.py’ was used cylindrical mesh of 10 µm length to simulate calcium bursts shown in Figure 5C.
6.	The script ‘HybridCaburst_stochSK.py’ was used with cylindrical mesh of 10 µm length to simulate calcium bursts shown in Figure 5D.
7.	The script ‘StochasticCaburst_cluster.py’ was used with cylindrical mesh of 40 µm length and variable number of CaP channels clustered with BK channels to simulate calcium bursts shown in Figure 6.
8.	The script ‘HybridCaburst_detchannels.py’ was used with all cylindrical meshes mentioned above to simulate calcium bursts and their analysis shown in Figure 7.
9.	The script ‘StochasticCaburst_wellmixed.py’ was used with all cylindrical meshes described above to simulate stochastic calcium bursts shown in Figure 8A-C. The races generated using this script were used for analysis in Figure 8D-F
10.	The script ‘StochasticCaburst_dendrite.py’ was used with dendritic mesh (‘chop_mesh.inp’) to simulate synaptically evoked dendritic calcium bursts shown in Figure 10C and Figure 11.
11.	The script ‘StochasticCaburst_dendrite_ampa.py’ was used with dendritic mesh (‘chop_mesh.inp’) to simulate dendritic calcium bursts shown in Figure 10D and Figure 12.

