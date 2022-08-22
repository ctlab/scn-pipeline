mkdir -p /mnt/space/conda-envs/
for YAML_FILE in /mnt/space/work/scn-pipeline/workflow/envs/*
do
echo $YAML_FILE
ENV=`basename $YAML_FILE .yaml`
mamba env create --prefix /mnt/space/conda-envs/$ENV --file $YAML_FILE
done