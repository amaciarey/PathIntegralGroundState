#!/bin/bash

# Reading some input parameters

PAR=$1
DEN=$2
POL=$3

# Check for correct number of parameters

if [ $# -ne 3 ] ; then
   echo 'The number of parameters is not correct!!!'
   echo ' '
   echo 'You must give the following input parameters:'
   echo '    1. Number of particles'
   echo '    2. Density of the system'
   echo '    3. Polarization angle'
   exit 0
fi

# Definition of some useful locations

BATCH='TRUE'
SRCDIR=$(pwd)'/PathIntegralGroundState'
CALDIR=$(pwd)'/calc'

WORKDIR=$CALDIR'/Np'$PAR'_n'$DEN'_a'$POL

# Create the working directory

if [ -e $WORKDIR ] ; then
   if [ -d $WORKDIR ] ; then
      cp $SRCDIR'/vpi' $WORKDIR
   else
      rm -f $WORKDIR 
      mkdir $WORKDIR
      cp $SRCDIR'/vpi' $WORKDIR
   fi
else
   mkdir $WORKDIR
   cp $SRCDIR'/vpi' $WORKDIR
fi

cd $WORKDIR

# Create the input file for the program

rm -f vpi.in
touch vpi.in

cat >vpi.in <<EOF
# RESTART PARAMETER                                                                     

.false.         # .true. Resume a previous simulation, .false. otherwise                
                                                                                       
# PHYSICAL DESCRIPTION OF THE SYSTEM                                                   
                                                                                       
2               # Dimensions                                                           
$PAR             # Number of particles                                                  
$DEN             # Density                                                              
$POL             # Polarization angle                                                   
.false.         # .true. if solid, .false. otherwise                                   
                                                                                       
# PARAMETERS THAT DESCRIBE THE MONTE CARLO SAMPLING                                     
                                                                                      
.false.         # .true. if diagonal simulation, .false. otherwise                    
3.d-5           # Time step                                                             
60              # Number of beads                                                       
1982            # Seed of the random number generator                                 
0.d0            # Parameters of the propagator decomposition                           
0.d0                                                                                  
0.2d0           # Size of the CM movements (in units of density^(1/dim))              
18              # Length of the Staging movements                                     
5               # Number of staging movements per step                                
                                                                                      
# PARAMETERS AFFECTING THE MODEL TRIAL WAVE FUNCTION                                  
                                                                                      
10000           # Grid points for the tabulation of the wave function                 
0.07d0          # Variational parameter                                               
.false.         # .true. if tabulate the wave function, .false. otherwise             
                                                                                      
# PARAMETERS FOR THE SAMPLING OF THE OBDM                                             
                                                                                      
700.d0          # 1st parameter for the OBDM Guess for importance sampling, \approx <K>
0.05d0          # 2nd parameter for the OBDM Guess for importance sampling, \approx n_0
10              # Number of evaluations of OBDM per step                              
5               # Maximum angular momentum partial wave in OBDM                       
                                                                                      
# PARAMETERS OF THE SAMPLING PROCESS                                                  
                                                                                      
2000            # Number of blocks                                                    
50              # Steps per block                                                     
40              # Points in the grid for g(r) and OBDM                                
20              # K points for S(k) evaluation                                        
EOF

# Create the batch file to submit the job to the queue 

if [ $BATCH = 'TRUE' ] ; then

   rm -f launch.sh
   touch launch.sh

   cat >launch.sh <<-EOF
#!/bin/bash

cd $WORKDIR
./vpi > vpi.out
EOF
   chmod u+x launch.sh

#   qsub launch.sh

else 

   echo "Starting the execution..."
   ./vpi

fi

