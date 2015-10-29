#!/usr/bin/python
# supercell.py
#
# Usage: supercell.py nx ny nz py.in
#
# Jack Deslippe

def main(argv = None):
   if argv is None:
      argv = sys.argv
   argc = len(argv)
   if argc != 6:
      print "\n   Usage: %s nx ny nz ecut local/nonlocal\n" % argv[0]
      return 2
   nx = int(argv[1])
   ny = int(argv[2])
   nz = int(argv[3])
   ecut = int(argv[4])
   blocal = argv[5]

   latVec = [[0 for col in range(3)] for row in range(3)]
   atmPos = [[0 for col in range(4)] for row in range(6)]
   atmPosOut = []
   nAtoms = 6

# Print Beginning of file

   print "&control"
   print "prefix = 'titania'"
   print "calculation = 'scf'"
   print "restart_mode = 'from_scratch'"
   print "wf_collect = .false."
   print "disk_io = 'none'"
   print "tstress = .false."
   print "tprnfor = .false."
   print "outdir = './'"
   print "wfcdir = './'"
   print "pseudo_dir = './'"
   print "/"
   print "&system"
   print "ibrav = 0"
   print "celldm(1) = 8.7671"
   print "nat = "+str(6*nx*ny*nz)
   print "ntyp = 2"
   print "nbnd = "+str(24*nx*ny*nz)
   print "ecutwfc = "+str(ecut)
   if (blocal == "local"):
     print "/"
   elif (blocal == "nonlocal"):
     print "input_dft = 'pbe0'"
     print "nqx1 = 1"
     print "nqx2 = 1"
     print "nqx3 = 1"
     print "exxdiv_treatment = 'gygi-baldereschi'"
     print "/"
   else:
     print "Neither local or nonlocal specified. Dying"

# Print Next Section

   print "&electrons"
   print "startingwfc='atomic'"
   print "electron_maxstep = 1"
   print "conv_thr = 1.0d5"
   print "mixing_mode = 'plain'"
   print "mixing_beta = 0.7"
   print "mixing_ndim = 8"
   print "diagonalization = 'david'"
   print "diago_david_ndim = 4"
   print "diago_full_acc = .true."
   print "/"
   print "CELL_PARAMETERS"

# Initial positions are now hardcoded
#
#   try:
#   h = open(fni, 'r')
#   r = h.readlines()
#   h.close()
#   for n in range(0,3):
#      s = r[n]
#      latVec.append(s.split())
#      print (n, latVec[n][0])
#   for n in range(3,len(r)):
#      s = r[n]
#      if (len(s.split()) == 4):
#         nAtoms = nAtoms + 1
#         atmPos.append(s.split())
#         print (n, atmPos[n-3][0],atmPos[n-3][3])
#
#   for i in range(len(latVec)):
#      for j in range(len(latVec[i])):
#         latVec[i][j] = float(latVec[i][j])
#
#   for i in range(len(atmPos)):
#      for j in range(1,len(atmPos[i])):
#         atmPos[i][j] = float(atmPos[i][j])

   latVec[0][0] = 1.0
   latVec[0][1] = 0.0
   latVec[0][2] = 0.0
   latVec[1][0] = 0.0
   latVec[1][1] = 1.0
   latVec[1][2] = 0.0
   latVec[2][0] = 0.0
   latVec[2][1] = 0.0
   latVec[2][2] = 0.638588

   atmPos[0][0] = 'Ti'
   atmPos[0][1] = 0.0
   atmPos[0][2] = 0.0
   atmPos[0][3] = 0.0
   atmPos[1][0] = 'Ti'
   atmPos[1][1] = 0.5
   atmPos[1][2] = 0.5
   atmPos[1][3] = 0.5
   atmPos[2][0] = 'O'
   atmPos[2][1] = 0.3050877
   atmPos[2][2] = 0.3050877
   atmPos[2][3] = 0.0
   atmPos[3][0] = 'O'
   atmPos[3][1] = -0.3050877
   atmPos[3][2] = -0.3050877
   atmPos[3][3] = 0.0
   atmPos[4][0] = 'O'
   atmPos[4][1] = 0.8050877 
   atmPos[4][2] = 0.1949123
   atmPos[4][3] = 0.5
   atmPos[5][0] = 'O'
   atmPos[5][1] = 0.1949123
   atmPos[5][2] = 0.8050877
   atmPos[5][3] = 0.5

# Debugging
#   print "\nNumber of atoms in input file: %i \n" % nAtoms

   nAtomsOut = 0
   atmPosTemp = [0,0,0]

   print nx*latVec[0][0], nx*latVec[0][1], nx*latVec[0][2]
   print ny*latVec[1][0], ny*latVec[1][1], ny*latVec[1][2]
   print nz*latVec[2][0], nz*latVec[2][1], nz*latVec[2][2]

   print "ATOMIC_SPECIES"
   print "Ti 47.867 Ti.pbe.nml"
   print "O 15.9994 O.pbe.nml"
   print "ATOMIC_POSITIONS crystal"

   for ix in range(0,nx):
      for iy in range(0,ny):
         for iz in xrange(0,nz):
            for n in range(0,nAtoms):
               nAtomsOut += 1
               atmPosTemp[0] = (atmPos[n][1] + ix)/nx
               atmPosTemp[1] = (atmPos[n][2] + iy)/ny
               atmPosTemp[2] = (atmPos[n][3] + iz)/nz
               atmPosOut.append([atmPos[n][0],atmPosTemp[0],atmPosTemp[1],atmPosTemp[2]])
               print atmPos[n][0],atmPosTemp[0],atmPosTemp[1],atmPosTemp[2]

   print "K_POINTS automatic"
   print "1 1 1 1 1 1"

if __name__ == "__main__":
   import sys
   import math
   import string
   sys.exit(main())
