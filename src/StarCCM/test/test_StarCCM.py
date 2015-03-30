
import unittest
import os
from StarCCM.StarCCM import StarCCM_wrapper


class StarCCM_wrapperTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    # add some tests here...

    def test_StarCCM_wrapper(self):
        sim = StarCCM_wrapper()

        sim.numNodes         =  5              #number of nodes to request
        sim.numCPUs          =  24             #number of CPUs per node to request
        sim.mpiProcs         =  24             #number of mpi processes
        sim.WallTime         = '2:00:00'       #wait time, format: 00:00:00
        sim.GroupID          = 'a0000'         #my group ID
        sim.jobName          = 'viper_starccm' #job name
        sim.queue            = 'devel'         #queue level

        sim.runVolGrid       = 1                 #generate volume grid. default to no(0)
        sim.javaBatch1File   = 'ATR_mesh'        #name of batch file for a meshing run without file extension .java
        sim.STLFile          = 'combined.stl'    #name of the STL file
        sim.simMeshFile      = 'g1_ATR_mesh.sim' #name of mesh file (a .sim file).

        sim.mySphereRadius   = 1500.0  #radius of freestream outer boundary in feet
        sim.mySphereX        = 0.0     #Center of freestream sphere in feet
        sim.mySphereY        = 0.0     #Center of freestream sphere in feet
        sim.mySphereZ        = 0.0     #Center of freestream sphere in feet
        sim.mySphereTriangles= 200.0   #size of sphere outer boundary facets
        sim.myPrismFieldRatio= 10.0    #thickness ratio of near field to outermost prism
        sim.myBLcells        = 9       #number of cells in boundary layer normal direction
        sim.myBLthickness    = 10.0    #thickness of boundary layer in inches
        sim.myBaseSize       = 8.0     #mesh base size in feet
        sim.myCurvature      = 90      #number of points to divide a circle
        sim.mySurfaceGrowthRate = 1.25 #growth rate (max size ratio) of surface triangles
        sim.myFeatureAngle   = 17.0    #maximum angle for defining sharp edges on ATR models in degrees
        sim.myMinMesh        = 0.2     # smallest cell size in percent
        sim.myEdgeTarget     = 0.4     # target size for feature curve edges in percent

        sim.makeSurfaceMesh  = True #use true to make surface mesh, false to skip
        sim.makeVolumeMesh   = True #use true to make volume mesh, false to skip
        sim.saveMeshFile     = True #use true to save final mesh file, false to skip

        sim.runCFD           = 1                 #run CFD. default to yes(1)
        sim.javaBatch2File   = 'ATR_flow'        #name of batch file for a CFD run without file extension .java
        sim.simFlowFile      = 'g1_ATR_flow.sim' #name of simulation file (a .sim file).

        sim.myMach           = 0.5     #Mach number
        sim.myAoA            = 3.0     #angle of attack, degree
        sim.myBeta           = 3.0     #beta, degree
        sim.myPressure       = 6.75886 #pressure, psi
        sim.myTemperature    = 447.415 #temperature, Rankine
        sim.myMoment_X       = 0.0     #moment reference location, X, feet
        sim.myMoment_Y       = 0.0     #moment reference location, Y, feet
        sim.myMoment_Z       = 0.0     #moment reference location, Z, feet

        sim.myCFL            = 5.0     #final CFL number
        sim.myIteration      = 5       #number of flow solver iterations
        sim.myCFLRampStart   = 0.1     #starting CFL number
        sim.myCFLRampIters   = 20      #number of iterations over which to ramp CFL up to final value
        sim.useExpert        = True    #flag whether use expert driver for solver control
        sim.myColorLevels    = 64      #number of colorbar levels in key

        sim.saveSimFile      = True #use true to save final sim file, false to skip
        sim.saveResidualImage= True #use true to save residuals image, false to skip
        sim.saveFandMImages  = True #use true to save final F&M images, false to skip
        sim.saveCpImages     = True #use true to save final Cp images, false to skip
        sim.saveResidualCSV  = True #use true to save residuals history, false to skip
        sim.saveFandMCSV     = True #use true to save final F&M histories, false to skip

        sim.cshUserBatch1    = 0                           #User will provide a shell script for mesh sim. default to no(0)
        sim.cshBatch1File    = 'startQsub_mesh_javaScript' #name of batch file for mesh sim without file extension .csh
        sim.cshUserBatch2    = 0                           #User will provide a shell script for CFD sim. default to no(0)
        sim.cshBatch2File    = 'startQsub_flow_javaScript' #name of batch file for CFD sim without file extension .csh
        sim.LicLocalPort     = '2999'                      #the port number on local machine to access license
        sim.LicAccessName    = 'bridge1'


        try:
            sim.run()
        finally:
            print "done"


if __name__ == "__main__":
    unittest.main()


