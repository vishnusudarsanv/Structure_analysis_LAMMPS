# Structure_analysis_LAMMPS
  This program analyzes the structure of polymer systems from LAMMPS trajectory
   files. It computes several structural properties including:
     - Solvation shell statistics for Li ions
     - End-to-end distance distributions for polymer chains
     - Persistence length via bond vector correlations
     - Radius of gyration (Rg) distributions
     - Radial distribution function (RDF) for polymer atoms

   The program reads atomic coordinates from a LAMMPS trajectory file, performs
   coordinate unwrapping for periodic boundary conditions (both orthorhombic and
   non-orthorhombic), and accumulates histograms and averages for the above
   properties over multiple simulation steps.

 Input Files:
   - LAMMPS trajectory file (wrap_mod.lammpstrj)

Output Files:
   - rg.dat: Histogram of radius of gyration
   - endtoend_4b_4.dat: Histogram of end-to-end distances
   - totrdf.dat: Radial distribution function data
   - persistence.dat: Persistence length data

 Key Parameters:
   - mxstep: Maximum number of simulation steps to analyze
   - nbins, rgbin, endnbins: Number of bins for RDF, Rg, and end-to-end histograms
   - natoms_peo: Number of atoms per PEO chain
   - totresid: Total number of polymer chains (residues)
   - nlip: Number of Li ions
   - Various bin sizes and geometric parameters for histogramming and box dimensions

 Major Variables:
   - coord_wrap, coord_unwrap: Wrapped and unwrapped coordinates for all atoms
   - li_coord: Coordinates of Li ions
   - hist_rdf_all, hist_rg, hist_enddist: Histograms for RDF, Rg, and end-to-end distance
   - t_bond, t_sep: Bond vectors and their correlations for persistence length
   - x_com, y_com, z_com: Center of mass coordinates for each polymer chain

 Usage Notes:
   - The program is tailored for a specific system (PEO/Li composite) and expects
     certain atom type conventions in the input trajectory.
   - Both orthorhombic and non-orthorhombic periodic boundary conditions are handled.
   - Output files are overwritten at each run.

!------------------------------------------------------------------------------
