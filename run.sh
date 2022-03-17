# input files
SDF_FILE="Enamine_Covalent_Screening_Library_11200cmpds_plated_20211006.sdf"
NAME="Enamine11000"

# prepare data
python src/process_sdf.py data/$SDF_FILE $NAME

# ersilia rdkit descriptors
ersilia serve rdkit-descriptors
ersilia -v api -i $NAME/smiles.smi -o $NAME/rdkit-descriptors.csv
ersilia close

# ersilia signaturizer
ersilia serve cc-signaturizer
ersilia -v api -i $NAME/smiles.smi -o $NAME/cc-signaturizer.h5
ersilia close

# ersilia grover
ersilia serve grover-embedding
ersilia -v api -i $NAME/smiles.smi -o $NAME/grover-embedding.h5
ersilia close

# diversity
python src/diversity.py $NAME

# report
python src/report.py $NAME