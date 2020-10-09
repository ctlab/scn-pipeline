import copy
import os
import pathlib
import re
import warnings

from jinja2 import Environment, FileSystemLoader

TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "templates")


def fill_template(j2_env, config: dict, file: str):
    return j2_env.get_template(file).render(**config)


def prepare_files(work_dir: str, config: dict) -> None:
    if not config['bulk_like']:
        with open(f"{work_dir}/Snakefile", "w") as snake_file, \
                open(f"{work_dir}/scripts/analysis.R", "w") as analysis_file, \
                open(f"{work_dir}/task.bash", "w") as task_file:
            j2_env = Environment(loader=FileSystemLoader(TEMPLATE_DIR),
                                 trim_blocks=True)
            task_file.write(fill_template(j2_env, config, 'task.bash'))
            snake_file.write(fill_template(j2_env, config, 'Snakefile'))
            analysis_file.write(fill_template(j2_env, config, 'analysis.R'))
        if config['AnalysisType'] == 'single':
            with open(f'{work_dir}/scripts/get_count_matrix.R', "w") as get_count_matrix_file:
                get_count_matrix_file.write(fill_template(j2_env, config, 'get_count_matrix.R'))
    else:
        with open(f"{work_dir}/Snakefile", "w") as snake_file, \
                open(f"{work_dir}/task.bash", "w") as task_file:
            j2_env = Environment(loader=FileSystemLoader(TEMPLATE_DIR),
                                 trim_blocks=True)
            task_file.write(fill_template(j2_env, config, 'task.bash'))
            snake_file.write(fill_template(j2_env, config, 'fluidigm_pipeline.smk'))



def prepare_rules(work_dir: str, config: dict) -> None:
    if not config['bulk_like']:
        j2_env = Environment(loader=FileSystemLoader(TEMPLATE_DIR),
                             trim_blocks=True)
        if config['AnalysisType'] == 'single':
            with open(f'{work_dir}/rules/get_data.smk', 'w') as get_f, \
                    open(f'{work_dir}/rules/run_kallisto.smk', 'w') as kallisto, \
                    open(f'{work_dir}/rules/get_counts.smk', 'w') as get_counts:
                get_f.write(fill_template(j2_env, config, 'rules/get_data.smk'))
                kallisto.write(fill_template(j2_env, config, 'rules/run_kallisto.smk'))
                get_counts.write(fill_template(j2_env, config, 'rules/get_counts.smk'))
            if config['cell_ranger']:
                with open(f'{work_dir}/rules/define_version.smk', 'w') as def_rule:
                    def_rule.write(fill_template(j2_env, config, 'rules/define_version.smk'))
                with open(f'{work_dir}/scripts/define_version.py', 'w') as dev_script:
                    dev_script.write(fill_template(j2_env, config, 'scripts/define_version.py'))
            else:
                with open(f'{work_dir}/scripts/prepare_kallisto.py', 'w') as prep_script:
                    prep_script.write(fill_template(j2_env, config, 'scripts/prepare_kallisto.py'))
        with open(f'{work_dir}/rules/secondary_analysis.smk', 'w') as seurat, \
                open(f'{work_dir}/rules/add_uns.smk', 'w') as uns, \
                open(f'{work_dir}/scripts/add_uns.py', 'w') as uns_py:
            seurat.write(fill_template(j2_env, config, 'rules/secondary_analysis.smk'))
            uns.write(fill_template(j2_env, config, 'rules/add_uns.smk'))
            uns_py.write(fill_template(j2_env, config, 'scripts/add_uns.py'))
    else:
        with open(f"{work_dir}/scripts/kallisto.sh", "w") as kallisto_script, \
                open(f"{work_dir}/scripts/srr_to_srs.R", "w") as srr_to_srs_script, \
                open(f"{work_dir}/scripts/srs_to_gse.R", "w") as srs_to_gse_script, \
                open(f"{work_dir}/rules/get_data.smk", "w") as get_data, \
                open(f"{work_dir}/rules/kallisto.smk", "w") as kallisto_rule, \
                open(f"{work_dir}/rules/srr_to_srs.smk", "w") as srr_to_srs_rule, \
                open(f"{work_dir}/rules/srs_to_gse.smk", "w") as srs_to_gse_rule:
            j2_env = Environment(loader=FileSystemLoader(TEMPLATE_DIR),
                                 trim_blocks=True)
            kallisto_script.write(fill_template(j2_env, config, 'scripts/fluidigm/kallisto.sh'))
            srr_to_srs_script.write(fill_template(j2_env, config, 'scripts/fluidigm/srr_to_srs.R'))
            srs_to_gse_script.write(fill_template(j2_env, config, 'scripts/fluidigm/srs_to_gse.R'))
            get_data.write(fill_template(j2_env, config, 'rules/fluidigm/get_data.smk'))
            kallisto_rule.write(fill_template(j2_env, config, 'rules/fluidigm/run_kallisto.smk'))
            srr_to_srs_rule.write(fill_template(j2_env, config, 'rules/fluidigm/srr_to_srs.smk'))
            srs_to_gse_rule.write(fill_template(j2_env, config, 'rules/fluidigm/srs_to_gse.smk'))



def analysis_run(config: dict) -> None:
    work_dir = config["AnalysisFolder"]
    pathlib.Path(work_dir).mkdir(parents=True, exist_ok=True)
    pathlib.Path(f'{work_dir}/rules').mkdir(parents=True, exist_ok=True)
    pathlib.Path(f'{work_dir}/scripts').mkdir(parents=True, exist_ok=True)
    prepare_files(work_dir, config)
    prepare_rules(work_dir, config)


def analysis(config: dict) -> None:
    pathlib.Path(config["PathToAnalysis"]).mkdir(parents=True, exist_ok=True)
    if config['Merge'] and not config['bulk_like']:
        assert len(set(config[
                           'Organism'])) == 1, "Experiment with samples from different organism: can't implement merge mode, try to separate yaml config by sample"
        if len(config['SampleIds']) > 1:
            config_local = copy.deepcopy(config)
            config_local["AnalysisFolder"] = os.path.join(config["PathToAnalysis"], "merged")
            config_local["AnalysisType"] = "many"
            config_local['FirstSample'] = re.sub('/counts.RData', '', config_local["Objects"][0])
            if not config_local["FilterUMI"] or config_local['WholeUMI']:
                analysis_run(config_local)
                config_local['WholeUMI'] = False
            if config_local["FilterUMI"]:
                config_local["AnalysisFolder"] = os.path.join(config_local["AnalysisFolder"], "filtered")
                analysis_run(config_local)
        else:
            warnings.warn("Less than two samples but `Merge` set True, skipping merge")

    if config['RunSingleAnalysis'] and not config['bulk_like']:
        for i, sample in enumerate(config["SampleIds"]):
            config_local = copy.deepcopy(config)
            config_local["AnalysisFolder"] = os.path.join(config["PathToAnalysis"],
                                                          config["SampleIds"][i])
            if config["Objects"]:
                config_local["Object"] = config["Objects"][i]
            else:
                config_local["PathToCount"] = config["PathsToCounts"][i]
            config_local["SampleId"] = config["SampleIds"][i]
            config_local["RunName"] = config["SampleIds"][i]
            config_local["AnalysisType"] = "single"
            config_local["Organism"] = config["Organism"][i]
            if not config_local["FilterUMI"] or config_local['WholeUMI']:
                analysis_run(config_local)
                config_local['WholeUMI'] = False
            if config_local["FilterUMI"]:
                config_local["AnalysisFolder"] = os.path.join(config_local["AnalysisFolder"], "filtered")
                analysis_run(config_local)

    if config['RunSingleAnalysis'] and config['bulk_like']:
        config_local = copy.deepcopy(config)
        config_local["AnalysisFolder"] = config["PathToAnalysis"]
        config_local["Organism"] = config["Organism"]
        analysis_run(config_local)