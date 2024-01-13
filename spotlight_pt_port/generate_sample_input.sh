#!/bin/bash
#
outf=$1
if test -z "$outf"
then
	outf="sample_input.dnvin"
fi

loc=`pwd`
cd ../dnv/data
datapath=`pwd`
cd ../../spotlight_pt_port

echo "Protein=path/to/protein.gro" > $outf
echo "ProteinFormat=gro" >> $outf
echo "Residues=1,3,4,17,35" >> $outf
echo "Seeds=All" >> $outf
echo "Dielectric=30" >> $outf
echo "RestrainDistance=0.6" >> $outf
echo "ProgramName=spotlight_jobname" >> $outf
echo "Sizes=16-26" >> $outf
echo "Step=1" >> $outf
echo "MolCount=100" >> $outf
echo "Oscillations=6" >> $outf
echo "MaxTries=500" >> $outf
echo "TriesPerLigand=75" >> $outf
echo "SeedCount=720" >> $outf
echo "Strategy_Restrain=Yes" >> $outf
echo "SourceFolder=$datapath" >> $outf
echo "Charges=FF" >> $outf
echo "Type=Standard" >> $outf
echo "Parallel=On" >> $outf
echo "Optimize=Yes" >> $outf

echo "Generated sample input file. Protein GRO file and residue numbers need to be entered"
