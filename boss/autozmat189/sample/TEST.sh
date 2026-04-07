#!/bin/sh

install_result () {
	echo " /bin/cp $scratch_file $AUTHENTIC_DIR"
	/bin/cp $scratch_file $AUTHENTIC_DIR
	if [ $? = 0 ]; then
		nsuccess=`expr $nsuccess + 1`
	else
		nfailure=`expr $nfailure + 1`
	fi

	echo ""
}

compare_result () {
	echo " diff $scratch_file $authentic_file > /dev/null 2>&1"
	diff $scratch_file $authentic_file > /dev/null 2>&1
	if [ $? = 0 ]; then
		echo " SUCCESS: $scratch_file and $authentic_file are identical"
		nsuccess=`expr $nsuccess + 1`
	else
		echo " ERROR: $scratch_file and $authentic_file are different"
		nfailure=`expr $nfailure + 1`
	fi
	echo ""
}

run_autozmat () {
	base=`echo $ifile | cut -f1 -d'.'`
	ofile=${base}.${ofmt}
	scratch_file=$SCRATCH_DIR/$ofile
	authentic_file=$AUTHENTIC_DIR/$ofile

	njobs=`expr $njobs + 1`
	echo "${njobs}. $message"
	echo " $AUTOZMAT -i $ifmt -o $ofmt $options $ATOMTYPE_OPT -f $scratch_file $dir/$ifile"
	$AUTOZMAT -i $ifmt -o $ofmt $options $ATOMTYPE_OPT -f $scratch_file $dir/$ifile > $errlog 2>&1

	if [ $? != 0 ]; then
		echo " ERORR: `basename $AUTOZMAT` failed."
		cat $errlog
		nfailure=`expr $nfailure + 1`
		return
	fi

	if [ "$INSTALL" = "TRUE" ]; then
		install_result
	else
		compare_result
	fi
}

#----------- MAIN -----------------------
AUTOZMAT=../bin/autozmat
AUTHENTIC_DIR=Authentic
SCRATCH_DIR=Scratch
ATOMTYPE_OPT=
INSTALL=FALSE

for arg in $*
do
	var=`echo $arg |cut -f1 -d'='`
	val=`echo $arg |cut -f2 -d'='`
	if [ "$var" = "" ]; then
		continue
	fi
	eval $var=$val
done

if [ "$ATOMTYPE_OPT" = "" ]; then
	AUTHENTIC_DIR=$AUTHENTIC_DIR/atomtyping
else
	AUTHENTIC_DIR=$AUTHENTIC_DIR/noatomtyping
fi

errlog=$SCRATCH_DIR/err.log
njobs=0
nsuccess=0
nfailure=0

#------------- BOSS Z-matrix -----------------------
dir=boss
ifmt=boss

ofmt=boss
ifile=acetaminophen.zmat
options=""
message="BOSS -> BOSS Z-Matrix
   with auto atomtyping"
run_autozmat

ofmt=mopac
ifile=acetaminophen.zmat
options=""
message="BOSS -> MOPAC Z-Matrix"
run_autozmat

#ofmt=pdb
#ifile=crotonate.zmat
#options="-p $dir/crotonate.par -x"
#message="BOSS Z-Matrix -> PDB
#   Dummies are removed and connection table is rebuilt."
#run_autozmat

ofmt=boss
ifile=crown.zmat
options="-z default2"
message="BOSS -> BOSS Z-matrix
   Auto atomtyping. Z-matrix options (non-explicit variable dihedrals)"
run_autozmat

ofmt=mmod
ifile=crown.zmat
options="-x"
message="BOSS Z-matrix -> MacroModel
   with connection table rebuilding"
run_autozmat

#------------- GAUSSIAN Z-matrix -----------------------
dir=gauss
ifmt=gauss

ofmt=boss
ifile=phacamide.com
options="-z default"
message="GAUSSIAN -> BOSS Z-Matrix
   with auto atomtyping and Z-matrix options"
run_autozmat

ofmt=jaguar
ifile=phacamide.com
options=""
message="GAUSSIAN -> JAGUAR Z-Matrix"
run_autozmat

ofmt=mmod
ifile=phacamide.com
options=""
message="GAUSSIAN Z-Matrix -> MacroModel"
run_autozmat

#------------- JAGUAR Z-matrix -----------------------
dir=jaguar
ifmt=jaguar

ofmt=gauss
ifile=etoh.in
options=""
message="JAGUAR -> GAUSSIAN Z-matrix"
run_autozmat

#------------- MCPRO Z-matrix -----------------------
dir=mcpro
ifmt=mcpro

ofmt=boss
ifile=diphe.zmat
options="-p $dir/diphe.par -z default -a"
message="MCPRO -> BOSS Z-matrix
   No atomtyping (-a) because atom types were already assigned in the Z-matrix."
run_autozmat

ofmt=pdb
ifile=diphe.zmat
options="-p $dir/diphe.par -x"
message="MCPRO Z-matrix -> PDB
   with dummies removed (-x)"
run_autozmat

ofmt=mdl
ifile=diphe.zmat
options="-p $dir/diphe.par -x"
message="MCPRO Z-matrix -> MDL MOL
   with dummies removed (-x)"
run_autozmat

ofmt=mol2
ifile=diphe.zmat
options="-p $dir/diphe.par -x"
message="MCPRO Z-matrix -> Tripos MOL2
   with dummies removed (-x)"
run_autozmat

ofmt=mind
ifile=diphe.zmat
options="-p $dir/diphe.par -x"
message="MCPRO Z-matrix -> MindTool
   with dummies removed (-x)"
run_autozmat

ofmt=mmod
ifile=diphe.zmat
options="-p $dir/diphe.par -x"
message="MCPRO Z-matrix -> MacroModel
   Dummies are removed (-x). Bond orders and connectivity are rebuilt."
run_autozmat

ofmt=pepz
ifile=diphe.zmat
options="-p $dir/diphe.par -a"
message="MCPRO Z-matrix -> PEPZ
   No atomtyping (-a) because atom types were already assigned in the Z-matrix."
run_autozmat

#------------- MDL -----------------------
dir=mdl
ifmt=mdl

ofmt=boss
ifile=chlorothiazide.mdl
options="-z default"
message="MDL -> BOSS Z-matrix
   Atomtyping and default Z-matrix options (-z default)"
run_autozmat

ofmt=mol2
ifile=chlorothiazide.mdl
options=""
message="MDL -> MOL2"
run_autozmat

ofmt=xyz
ifile=chlorothiazide.mdl
options=""
message="MDL -> XYZ"
run_autozmat

#------------- MindTool -----------------------
dir=mind
ifmt=mind

ofmt=mmod
ifile=catolf.mind
options="-a -g break_mind"
message="MindTool -> MacroModel
   Break '%' separated atoms into individual molecules (-g break_mind)"
run_autozmat

ofmt=pdb
ifile=3cro.mind
options="-a"
message="MindTool -> PDB"
run_autozmat

#------------- MacroModel -----------------------
dir=mmodel
ifmt=mmod

ofmt=mdl
ifile=cyclododecane.dat
options=""
message="MacroModel (multiple structures) -> MDL SD"
run_autozmat

ofmt=mind
ifile=dynamics.dat
options=""
message="MacroModel (dynamics) -> MindTool"
run_autozmat

#------------- MOL -----------------------
dir=mol
ifmt=mol

ofmt=mdl
ifile=
options=""
message="Tripos MOL -> MDL"
#run_autozmat

#------------- MOL2 -----------------------
dir=mol2
ifmt=mol2

ofmt=mdl
ifile=chlorothiazide.mol2
options=""
message="Tripos MOL2 -> MDL"
run_autozmat

ofmt=boss
ifile=dephenhydramine.mol2
options="-z default"
message="Tripos MOL2 -> BOSS Z-matrix
   Auto atomtyping and default Z-matrix options"
run_autozmat

ofmt=mcpro
ifile=viagra.mol2
options="-z default"
message="Tripos MOL2 -> MCPRO Z-matrix
   Auto atomtyping and default Z-matrix options"
run_autozmat

#------------- MOPAC -----------------------
dir=mopac
ifmt=mopac

ofmt=gauss
ifile=tetrazine.dat
options=""
message="MOPAC -> GAUSSIAN Z-Matrix"
run_autozmat

#------------- PDB -----------------------
dir=pdb
ifmt=pdb

ofmt=mmod
ifile=3cro.pdb
options="-a"
message="PDB -> MacroModel"
run_autozmat

ofmt=mcpro
ifile=barnase.pdb
options="-a -z default"
message="PDB -> MCPRO Z-matrix
   Atomtyping and default Z-matrix options"
run_autozmat

#------------- XYZ -----------------------
dir=xyz
ifmt=xyz

ofmt=boss
ifile=ph2acetyl.xyz
options="-z default"
message="XYZ -> BOSS Z-matrix
   Auto atomtyping and default Z-matrix options"
run_autozmat

ofmt=xyz
ifile=benzene.xyz
options="-H"
message="XYZ (benzene with no hydrogens) -> PDB (benzene with hydrogens)"
run_autozmat

ofmt=pepz
ifile=ph2acetyl.xyz
options="-g pepz=500"
message="XYZ -> PEPZ
   Atom types start from 500."
run_autozmat

#######################################################
###   SUMMARY
#######################################################

success_percent=0
failure_percent=0
if [ "$njobs" -gt 0 ]; then
	success_percent=`expr $nsuccess \* 100 / $njobs`
	failure_percent=`expr $nfailure \* 100 / $njobs`
fi

echo "---------------------------------------------------------------"
echo "                SUMMARY OF TEST RESULTS                        "
echo "---------------------------------------------------------------"
echo " Number of Succeeded Tests        =     $nsuccess ($success_percent%)"
echo " Number of Failed Tests           =     $nfailure ($failure_percent%)"
echo "---------------------------------------------------------------"
echo " Toal Number of Tests             =     $njobs"

exit $nfailure

