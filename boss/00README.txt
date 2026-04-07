
 
                     BBBBB     OOOO     SSSS     SSSS
                     B    B   O    O   S    S   S    S
                     B    B   O    O   S        S
                     BBBBB    O    O    SSSS     SSSS
                     B    B   O    O        S        S
                     B    B   O    O   S    S   S    S
                     BBBBB     OOOO     SSSS     SSSS

                      W. L. Jorgensen, Yale University

                       Version 5.1, April 2024 Release

Some basic information to get started with BOSS under UNIX or linux:

    USER'S MANUAL
    The PDF file for the BOSS user's manual is in
    bossman/boss51.pdf
    It can be displayed with Adobe Acrobat Reader.                        
    Minimally, read page 6 of the user's manual (Can't Wait to Start).
    And, see the New Features section, page 9 in the manual, for a
    description of the important new capabilities in BOSS.

    For UNIX or Linux computers, 
    define the path to the top BOSS directory in your .cshrc
    file with an entry like
    setenv BOSSdir ~/boss

    For Linux, please note that all provided scripts that just execute BOSS once, 
    such as  xOPT, execute BOSS with a command like:
    csh commandfile >& log  
    For some test jobs that perform iterative job submissions, the user may prefer
    to submit a command file with the -e flag (exits if there is an error 
    without further processing):
    csh -e commandfile >& log &
    However, this will not be allowed on your system if the default system start-up
    file, /etc/csh.cshrc, causes an error status of any kind. The remedy is to not
    use the -e flag or comment out the offending lines in /etc/csh.cshrc
  
    If you are an administrator, you may wish to skip to the
    installation instructions below. 
-----------------------------------------------------------------------------------------
    Many scripts are provided for execution of common jobs like energy minimizations,
    conformational search, semiempirical MO calculations, and Monte Carlo simulations
    for a solute in water are provided in the scripts directory.
    See the 00README file there.

    Use of BOSS should be acknowledged as
    Jorgensen, W. L.; Tirado-Rives, J. J. Comput. Chem. 2005, 26, 1689-1700.

    Key New Features in BOSS 5.1:
    (1) Added LJ cutoff correction for solutes in any solvent; this is automatic. 
    (2) Added ability to write zmats for Charge and LJ annihilations. 
    (3) This is used in new test job FEPannihil for automatic calculation of free
        energies of hydration - see JPCB 2024, 128, 250-262.
    (4) Increased output formats for longer runs. Averaging can now be as long
        as 9,999,999,999 configurations. 
    (5) Test job MCsolute was added for automatic execution of MC simulations 
        for a solute in TIP4P water for 0.5 B or 5 B configurations. This allows
        direct calculations of heats and volumes of hydration.  
        See PCCP 2024, 26, 8141-8147.
    (6) Directory boss/molecules/OPLS2020 now contains many z-matrices for 
        molecules with OPLS/2020 type assignments.

    Key New Features in BOSS 5.0:
    (1) Implementation of 1.14*CM1A-LBCC charges for neutral molecules, and the 
        scripts xAM1SPL, xAM1OPTL, and xZCM1AL.
    (2) Addition of argon N=512 box to allow creation of custom solvents with N<=512
    (3) Addition of MixSolvent test job; example is water saturated octanol.
    (4) The OPLS2020 parameters for alkanes are provided. New z-matrices are in 
        molecules/small/zzAlkanes2020
    (5) Cation-pi interactions can be treated with 1/r^4 terms added to the force field.
        See JCTC 2020, 16, 7184-7194, and test job cationpi.  

    Key New Features in BOSS 4.9:
    (1) Treatment of halogen bonding has been included in the force field by addition
        of partial positive point charges on the C-X axis for X = Cl, Br, and I. See the
        new testjob HalogenBond for examples. New atom types XC, XB, and XI were needed.
    (2) New OPLS-AA alkane parameters have been provided - see types 54 - 61 in 
        oplsaa.par. These give improved liquid properties for branched alkanes.
        See molecules/small/zzNewAlkanes for examples of new Z-matrices.
    (3) A utility program typeQLJ has been provided and can be executed with the
        script xQLJ. The input is the PDB file for a molecule; output is the 
        BOSS atom types, atomic charges (OPLS/CM1A) and Lennard-Jones parameters.
    (4) Several new atom types (OA, SA, NS, NX, CF) were added - see notes/atomtypes.txt
    (5) Torsional energetics have been improved for biaryl and sustituted heteroaryl
        systems.
    (6) The OPLS-AA/M force field has been implemented for peptides and proteins. 

    Key New Features in BOSS 4.8:
    (1) Overlap sampling FEP calculations have been implemented; see the new SOSgas and
        SOSaq testjobs and the FEPreview preprint in the bossman directory.
    (2) The xFEPgas and xFEPaq test jobs replace the former fepme and have been expanded to 
        illustrate many PhX to PhY FEP calculations.
    (3) IRIX executables are no longer included.
    (4) Many additions have been made to the OPLS-AA force field; replace any old copies
        of oplsaa.par and oplsaa.sb with the new provided ones. The parameters for
        alkali and halide ions were optimized - see the Ions2006 reprint in bossman.
    (5) Treatment of conjugated systems especially trienes
        has improved with addition of atom type C# for C3 and C4 of a triene.

    Key New Features in BOSS 4.7:
    (1) Solute-solute and solvent by solute polarization are optionally treated for
        all calculations including energy minimizations and MC/FEP calculations.
        This is potentially particularly important for cation-pi interactions and
        polar solutes in low-dielectric media.
        See the new testjob, feppolrz, for details.
    (2) See the new scripts xPOPT, xPOPTF, and xZCM1A
    (3) Treatment of allenes has been added - see allen.z, allen1m.z, and
        allen2m.z in molecules/small
    (4) Improvements were made for the optimizations of complexes using the
        xOPT and xOPTF scripts. Resubmission is no longer needed with either. 
    (5) Test job etoh has been added. It is an example of a MC simulation for
        a solute in a custom solvent - ethanol in liquid ethanol.

    Key New Features in BOSS 4.6:
    (1) Implementation of the GB/SA hydration model for OPLS-AA optimizations, single point
        calculations, dihedral driving, conformational search, and standard MC simulations.
        Try the xSPGB, xCSGB100, xMCGB, xOPTGB,
        xOPTFGB (recommended for optimizations) and other GB scripts.
    (2) Extension of PDDG/PM3 to sulfur, silicon & phosphorus: JCTC 2005, 1, 817-823.          
    (3) Atomic charges from quantum calculations (CM1, CM3, Mulliken) are symmetrized for XHn
        groups and mono- and di-substituted benzenes.
    (4) The default scaling factor for CM1A charges has been optimized to be 1.14. This is
        now applied when using xAM1SP and xAM1OPT scripts. This is the recommended general
        charge model. All Z-matrices in the drugs directory have been redone with the
        symmetrized 1.14*CM1A charges. Please see the CMcharges.pdf file in the bossman
        directory for a preprint.
    (5) See the user's manual and the 00README file in the scripts directory about new scripts.
    (6) The MUSIC test job has been added - for multiple copy simulations for a probe molecule
        normally with a protein. The maximum number of solute atoms has been increased to 3500
        and of residues to 395.
    (7) All molecule Z-matrices are now in the molecules directory.

    ********** Job Scripts ********
    Scripts for many standard jobs are provided in the scripts directory. This has
    been expanded. Please see the scripts/00README file and the user's manual, page 5.

    ****** Conformational Search *******
    Conformational search (CS) is one of BOSS' most popular features - see the  
    consearch test job for more details.

    ********** Molecular Structures **********
    Z-matrices are in the molecules directory.
    Many new Z-matrices have been added to the small and optional drugs,
    antib, and drugHIV subdirectories.

    ********** Larger Molecules *********
    Energy minimizations and conformational search can now be performed on molecules 
    with up to 450 degrees of freedom (variables). 
    AM1 and PM3 calculations can be performed on molecules with up to 100
    non-H and 150 H atoms. For molecules with >100 atoms, the output Z-matrix
    with CM1 charges now uses atom types 800-899, as before, and then adds the remainder
    with types 9500 and up, e.g.
...
 897  1 HC   0.082644  2.500000  0.030000
 898  6 CT  -0.175417  3.500000  0.066000
 899  1 HC   0.110032  2.500000  0.030000
9500  1 HC   0.066297  2.500000  0.030000
9501  1 HC   0.067625  2.500000  0.030000
9502  1 HC   0.083185  2.500000  0.030000
...
    BOSS >4.3 recognizes as input both this format and the old format that only has 3-digits 
    for the atom type numbers. The first entry needs to be type 800 (800 for the old format
    or 0800 for the new).
    
-----------------------------------------------------------------------------------------

                          General Notes and Installation

(1) The default executables are for Linux and were obtained with the Intel FORTRAN
    compiler on a 2.66 GHz Intel Core2i quad running Fedora core release 12.  
    The -static flag has been used to
    include all libraries in the executables. They should work on
    any Linux system with the 2.4 or newer kernel.

    Please also see the information on autozmat in point (12) below.

    *** IMPORTANT ***
(2) Again, you need to define the path to the top BOSS directory in your .cshrc
    file with an entry like
    setenv BOSSdir ~/boss 
    *****************

    You do not have to do the following, but
    for general use of BOSS, we place the executable BOSS* in the
    usr/local/bin directory and the stored solvent boxes (watbox,
    org1box, and org2box) in usr/local/lib, so individual users do
    not need their own copies of these files. 

    The provided cmd files presume BOSS, watbox, org1box, org2box, and the
    parameter files are in BOSSdir. You will need to modify the command files,
    if you place these files in other directories.

    All of the scripts like xOPT in the boss/scripts directory could
    also be placed in /usr/local/bin so they can be executed from
    any directory.

(3) Files for the many test jobs have been provided. It is strongly
    recommended that the user's manual be read, and the test jobs
    gone over carefully and executed.
  
    Start with test job 1 and perform some gas-phase energy minimizations.

    The files for the test jobs are in the subdirectories of testjobs.
    The test jobs cover most typical applications of BOSS.
    New applications can generally be modeled
    on a test job; it is typical to begin a new project by creating a
    new directory and copying all of the files from the related test job
    to the new directory and then modifying the zmat, par and cmd files.

(4) If a user just wishes to use the standard scripts for optimizations,
    conformational search, etc., just make a new directory containing the 
    input structure files, copy the appropriate scripts to it from
    the boss/scripts directory, and execute them in the new directory.
    Or, place the scripts in usr/local/bin so they can be executed from
    any directory.

    Otherwise, for a new project, a user just needs a command file (cmd), a
    Z-matrix file (zmat), a parameter file (par), and maybe a coordinate
    input file (in) in their directory. In the usual case for a new
    project, the in file is empty.

    We always obtain the cmd and zmat files by just copying
    those files from a related test job or prior project, and then
    renaming and editing them for the new project.
    A new par file is readily the same way or by using the
    interactive pargen program. Just enter pargen.

    The locations of the BOSS*, watbox, org1box, org2box, and .sb
    files need to be indicated correctly in the cmd file. They are
    assumed to be in the directory defined by the setenv BOSSdir
    command.

(5) When the four files are ready, the program may be executed by
    csh xxxcmd >& xxxlog &
    where xxx indicates the new project prefix.

(6) The key thing with BOSS is having correct Z-matrices. 
    The Z-matrix file is critical because
    it defines the structure & atom types and
    for an FEP calculation how they change from reference to
    perturbed structure, and what variables are to be sampled
    in an MC simulation or optimized in an energy minimization.
    BOSS does contain a routine that checks for common mistakes
    in the solute Z-matrix file; if potential problems are found,
    warnings are provided in the ot file - check for them when
    a new Z-matrix is first used. To have a solvent Z-matrix checked,
    run a gas-phase optimization on it as a solute. 

    Hundreds of all-atom Z-matrices can be found in the 
    molecules/small subdirectory.

    A Z-matrix for each amino acid is provided in subdirectory peptide.

    Z-matrices for many drugs are provided in subdirectory drugs.

    The autozmat utility can be used to generate good BOSS Z-matrices
    from many other file types.

(7) The provided OPLS parameter files are extensive for
    atom types and stretching & bending parameters. 
    Sometimes there are missing torsional parameters, which
    may require some fitting for new systems. However,
    the autotyping and matching procedures for dihedral angles in
    BOSS are powerful and are adequate for most systems.

    The oplsaa.par and oplsaa.sb files contain the latest and complete
    parameters for the OPLS-AA force field.  
    *********************************************************************
    The all-atom parameter files, oplsaa.par and oplsaa.sb, are
    privileged information since much of it has not been published.
    These are the latest parameters for the OPLS-AA force field. The
    oplsaa.par file is normally automatically appended to make up the bottom 
    section of a par file for BOSS; the top section can be created by pargen.
    A paper containing many OPLS-AA parameters is in
    J. Am. Chem. Soc. 118, 11225-11242 (1996). See oplsaa.par for
    more references.

    The parameters are not to be distributed (or is any other part of the
    BOSS programs and files) without permission from W. L. Jorgensen.

    The all-atom files provide coverage for a wide range of organic
    and biomolecular systems. 
    *********************************************************************

    The subdirectory solbox contains additional Z-matrix and in files
    for equilibrated boxes of custom solvents.

(8) The subdirectory notes contains some additional information on
    specific topics. In particular, aapar.note specifies the OPLS-AA
    atom types for amino acids.

(9) The executable autozmat provides for interconversion of many file
    types. For example, it can create a BOSS Z-matrix file from a
    Sybyl Mol2 file or an MDL mol file. 
    Enter "autozmat -h" if you have any questions. Standard usage for
    conversion of a Sybyl mol2 file to a fully flexible BOSS Z-matrix is:
    ~/boss/autozmat -i mol2 -z default < file.mol2 > file.z
    This is simplified with provided scripts, e.g.
    xMOL2Z zoloft
    converts the Tripos mol2 file, zoloft.mol2, to zoloft.z, and
    xMOLZ zoloft
    converts the MDL mol file, zoloft.mol, to zoloft.z. 
    zoloft.mol2 and zoloft.mol can be found in the 
    molecules/small subdirectory.
    autozmat was developed by Dr. Dongchul Lim.

    Please also check our website, http://jorgensenresearch.com
    for recent developments and publication lists.
