#!/bin/bash

# Matthew Lee
# Fall 2019
# JagResearch
# MultiMutant

# SCRIPT ARGS:
# pdbID
# chainID
# resNumRange (X:Y)
# numMutations

# FLAGS:
# -em (energy minimization)
# -hphilic (mutate to only hydrophilic amino acids)
# -hphobic (mutate to only hydrophobic amino acids)

# CHANGES TO MAKE:
# count args/warn for incorrect flags/args
# rename wl type to something in the vein of "WL_$1".pdb
# add flag for logging output to a file

###### VARIABLE AND FUNCTION DECLARATIONS/DEFINITIONS ######

### Amino Acids ###

# all amino acids
declare -a allAcids=("A" "R" "N" "D" "C" "Q" "E" "G" "H" "I" "L" "K"
                                    "M" "F" "P" "S" "T" "W" "Y" "V")
# hydrophilic amino acids
declare -a hPhilicAcids=("S" "T" "N" "Q" "C" "E" "R" "H" "D")
# hydrophobic amino acids - gly, pro, phe, ala, Ile, Leu, Val
declare -a hPhobicAcids=("G" "P" "F" "A" "I" "L" "V" )


### Flag Options ###

# which amino acids to mutate to
declare -a amAcids=("${allAcids[@]}")
# energy minimization argument passed to rMutant
# (defaults to no, change with -em flag)
ENERMIN="no"

### Functions ###

# runs rMutant with the given arguments
# ARGS: pdbID chainID resNum mutTarget (optional)energyMinimization
runMutant () {
    #need to find way to track what forks have died to perform file operations
    #after all rMutants have finished
    #For debugging add this to end of script: && echo "Mutation of $3 to $4 complete." &
    ./proMute $1 $2 $3 $4 $5 &>/dev/null &
        #& >/dev/null && mkdir ../$6_out/$1.$2$3$4 && mv $1.$2$3$4.pdb $1.$2$3$4_em.pdb $1.$2$3$4.fasta.txt ../$6_out/$1.$2$3$4 -f 2>/dev/null &
}


#calls runMutant for each amino acid type
# ARGS: pdbID chainID resNum em foldername
allTargets() {
    for j in "${amAcids[@]}"
    do
                # ARGS: pdbID chainID resNum newAmino em foldername
        runMutant $1 $2 $3 $j $4
    done
}


###### EXECUTED PART ######
startTime=$(date +%s)

# Download proMutant/rMutant if it doesn't exist
if [ ! -d promute ];
then
        git clone https://gitlab.cs.wwu.edu/carpend3/promute.git
        cd promute
        make
        cd ..
        mv promute promute-src
        mv promute-src/build promute
        rm -rf promute-src
        
fi

#handle flags
for flag in "$@"
do
        if [ $flag == "-em" ];
        then
                ENERMIN="em"
        elif [ $flag == "-hphilic" ]
        then
                amAcids=("${hPhilicAcids[@]}")
        elif [ $flag == "-hphobic" ]
        then
                amAcids=("${hPhobicAcids[@]}")
        fi
done

#download specified protein pdb file from rcsb
if [ ! -f "promute/${1}.pdb" ]
then
        wget -q -O "promute/${1}.pdb" "https://files.rcsb.org/download/${1}.pdb"
fi

#grabbing range of residues (inclusive)
beginning=${3%:*}
end=${3#*:}

#calculate a 0-indexed number of mutations
numMutations=`expr $4 - 1`

#make output folder if it doesn't exist
if [ ! -d $1$2$3_$4_out ];
then
        mkdir $1$2$3_$4_out
fi

#make output tmp directory if it doesn't exist
if [ ! -d $1$2$3_$4_out/tmp ];
then
        mkdir $1$2$3_$4_out/tmp
fi

#copying pipeline invocation script into output folder 
cp sequentialPipelineInvocation.sh $1$2$3_$4_out/sequentialPipelineInvocation.sh

#cd-ing into folder containing rMutant
cd promute

count=0
#loop for each round of mutations.
for (( i=0; i<=$numMutations; i++ ));
do      
        echo "Starting mutation round $i"
        FILES=*.pdb

        # Combinations of mutations are restricted to strictly increasing sequences
        # All combinations are captured, assuming order doesn't matter
        lowerLim=`expr $beginning + $i`
        upperLim=`expr $end + $i - $numMutations`
        for (( j2=$lowerLim; j2<=$upperLim; j2++));
        do
                for f in $FILES;
                        do

                        

                        #Check if this file is from the last round of mutations by counting the "." in its name
                        IFS='.' read -ra splitFileName <<< "$f"
                        if [ "${#splitFileName[@]}" -eq $((i+2)) ]
                        then

                                prevResNum=0
                                if [ $i -ne 0 ]
                                then 
                                        # If this is not the first round of mutations, get the most recenly mutated resNum
                                        #get most recently mutated residue from filename
                                        prevMutation=${splitFileName[-2]}
                                        prevResNum=${prevMutation:1:-1}
                                fi

                                # Only mutate residues in ascending order
                                if [ $prevResNum '-lt' $j2 ];
                                then
                                        if [ $i -eq $numMutations ]
                                        then
                                                # Last round of mutations, call with em flag
                                                allTargets ${f::-4} $2 $j2 $ENERMIN
                                                count=$((count+1))
                                        else
                                                # Not the last round of mutations, no em flag
                                                allTargets ${f::-4} $2 $j2 "no"
                                        fi
                                fi
                        fi
                done
        done
        wait
done


echo "Finished $count mutations , cleaning..."
#Move all fasta.txt and .pdb files in /promute to the output folder, move WT pdb back to /promute.
FILES=*"_em.pdb"
for f in $FILES
do
    base=${f::-7}
    GROUP=$base*
    mkdir ../$1$2$3_$4_out/$base&>/dev/null && mv $GROUP ../$1$2$3_$4_out/$base&>/dev/null &
done
wait

FILES=*.fasta.txt
for f in $FILES
do
        
        if [ "$f" != "*.fasta.txt" ]
        then
                base=${f::-10}
                GROUP=$base*
                mkdir ../$1$2$3_$4_out/tmp/$base&>/dev/null && mv $GROUP ../$1$2$3_$4_out/tmp/$base&>/dev/null &
        fi
done
wait

endTime=$(date +%s)
endTime=`expr $endTime - $startTime`
echo "Runtime: $endTime Seconds"