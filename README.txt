This is a wrapper of Star-CCM+. It can generate both volume grid and CFD solutions.
User is expected to have access to Star-CCM+ executable and power-on-demand license.
It is also expected that the code will be run in PBS. The following system variables
are expected:

LM_LICENSE_FILE (e.g. 1999@localhost)
STAR_POWER_ON_DEMAND_LIC (e.g. ADeged4oUsdf6ghTln7ALP)

If the user wants to generate a volume grid (runVolGrid = 1), an surface grid file must 
have been provided:

- STLFile (a .stl file)
 
By default, two more input files are automatically generated: 

- javaBatch1File: a Star-CCM+ macro (a .java file),
- cshBatch1File:  a C-shell script (a .csh file).

However, user can choose to provide his/her own C-shell script (cshUserBatch1 = 1).
If everything goes smoothly, it will generate a volume grid:

- simMeshFile: the volume grid (a .sim file).

Volume grid can be used to generate CFD solutions.

If the user wants to generate a CFD solution (runCFD = 1), a mesh file (the simMeshFile 
file) must have been either generated in the previous step or provided by the user.

By default, two more input files are automatically generated:

- javaBatch2File: a Star-CCM+ macro (a .java file),
- cshBatch2File:  a C-shell script (a .csh file).

However, user can choose to provide his/her own C-shell script (cshUserBatch2 = 1).
If everything goes smoothly, it will generate a CFD solution:

- simFlowFile: the CFD solution (a .sim file).

several CSV and PNG files can be generated too.

The inputs are in imperial system.

