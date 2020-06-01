#Calls all aggregate py scripts for all mutants

SCRIPTS=aggregate_py/*
FILES=preliminary_exp_output/additional_EM_logs/*

for f in $FILES;
do
    IFS='/' read -ra splitPathName <<< "$f"
    IFS='_' read -ra splitFileName <<< "${splitPathName[2]}"

    for s in $SCRIPTS;
    do
        python ./$s ${splitFileName[0]}
        python ./$s ${splitFileName[0]} -b
    done
done    

