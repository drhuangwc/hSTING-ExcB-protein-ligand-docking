# hSTING-ExcB-protein-ligand-docking


hSTING-ExcB protein-ligand docking


[![DOI](https://zenodo.org/badge/648071283.svg)](https://zenodo.org/badge/latestdoi/648071283)


The three-dimensional structures of ligands for the present study were generated using Python programming language and the open-source RDKit library. The protein structure of full-length hSTING for docking experiments was obtained from the Protein Data Bank (PDB) with PDB ID of 6NT5. The protein structure was prepared and protonated with AMBER-FB15 force field. Protein-ligand blind docking was conducted to analyze the potential binding pocket and the protein-ligand interactions using a three-step method. The cavities of protein surface were firstly determined with CB-Dock. Second, the AutoDock Vina was employed to dock ligands following the protein-ligand docking protocol. Finally, based on the previous blind docking result, a covalent docking on Cys91 residue was conducted with explicitly specified binding site flexibility for four arginines residues (Arg83, Arg86, Arg94, and Arg95) using the AutoDockFR with the AutoDock4 scoring function following the covalent docking protocol. The covalent docking poses of the hSTING-Cys91-ExcB protein-ligand structure with lowest energy were obtained from the analysis. To further calculate the membrane embedding structure of the hSTING-Cys91-ExcB, the Bilayer Builder function in the CHARMM-GUI was performed, where PPM 2.0 was utilized for protein-membrane orientation calculation, and heterogeneous lipids with simplified composition of mammalian endoplasmic reticulum membranes, were utilized to simulate the membrane structure.
