import copy
import os
import pathlib
import warnings

from jinja2 import Environment, FileSystemLoader

TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "templates")

def fill_template(j2_env, config: dict, file: str):
    return j2_env.get_template(file).render(**config)


def prepare_files(work_dir: str, config: dict) -> None:
    with open(f"{work_dir}/Snakefile", "w") as snake_file, \
            open(f"{work_dir}/scripts/get_count_matrix.R", "w") as get_count_matrix_file, \
            open(f"{work_dir}/scripts/analysis.R", "w") as analysis_file, \
            open(f"{work_dir}/task.bash", "w") as task_file:
        j2_env = Environment(loader=FileSystemLoader(TEMPLATE_DIR),
                             trim_blocks=True)
        task_file.write(fill_template(j2_env, config, 'task.bash'))
        snake_file.write(fill_template(j2_env, config, 'Snakefile'))
        get_count_matrix_file.write(fill_template(j2_env, config, 'get_count_matrix.R'))
        analysis_file.write(fill_template(j2_env, config, 'analysis.R'))


def prepare_rules(work_dir: str, config: dict) -> None:
    j2_env = Environment(loader=FileSystemLoader(TEMPLATE_DIR),
                         trim_blocks=True)
    if config['cell_ranger']:
        with open(f'{work_dir}/rules/define_version.smk', 'w') as def_rule:
            def_rule.write(fill_template(j2_env, config, 'rules/define_version.smk'))
        with open(f'{work_dir}/scripts/define_version.py', 'w') as dev_script:
            dev_script.write(fill_template(j2_env, config, 'scripts/define_version.py'))
    else:
        with open(f'{work_dir}/scripts/prepare_kallisto.py', 'w') as prep_script:
            prep_script.write(fill_template(j2_env, config, 'scripts/prepare_kallisto.py'))
    with open(f'{work_dir}/rules/get_data.smk', 'w') as get_f, \
            open(f'{work_dir}/rules/run_kallisto.smk', 'w') as kallisto, \
            open(f'{work_dir}/rules/get_counts.smk', 'w') as get_counts, \
            open(f'{work_dir}/rules/secondary_analysis.smk', 'w') as seurat, \
            open(f'{work_dir}/rules/add_uns.smk', 'w') as uns, \
            open(f'{work_dir}/scripts/add_uns.py', 'w') as uns_py:
        get_f.write(fill_template(j2_env, config, 'rules/get_data.smk'))
        kallisto.write(fill_template(j2_env, config, 'rules/run_kallisto.smk'))
        get_counts.write(fill_template(j2_env, config, 'rules/get_counts.smk'))
        seurat.write(fill_template(j2_env, config, 'rules/secondary_analysis.smk'))
        uns.write(fill_template(j2_env, config, 'rules/add_uns.smk'))
        uns_py.write(fill_template(j2_env, config, 'scripts/add_uns.py'))

def analysis_run(config: str) -> None:
    work_dir = config["AnalysisFolder"]
    pathlib.Path(work_dir).mkdir(parents=True, exist_ok=True)
    pathlib.Path(f'{work_dir}/rules').mkdir(parents=True, exist_ok=True)
    pathlib.Path(f'{work_dir}/scripts').mkdir(parents=True, exist_ok=True)
    prepare_files(work_dir, config)
    prepare_rules(work_dir, config)


def analysis(config: str) -> None:
    pathlib.Path(config["PathToAnalysis"]).mkdir(parents=True, exist_ok=True)
    if config['Merge']:
        assert len(set(config[
                           'Organism'])) == 1, "Experiment with samples from different organism: can't implement merge mode, try to separate yaml config by sample"
        if len(config['SampleIds']) > 1:
            config_local = copy.deepcopy(config)
            config_local["AnalysisFolder"] = os.path.join(config["PathToAnalysis"], "merged")
            config_local["AnalysisType"] = "many"
            if not config_local["FilterUMI"] or config_local['WholeUMI']:
                analysis_run(config_local)
                config_local['WholeUMI'] = False
            if config_local["FilterUMI"]:
                config_local["AnalysisFolder"] = os.path.join(config_local["AnalysisFolder"], "filtered")
                analysis_run(config_local)
        else:
            warnings.warn("Less than two samples but `Merge` set True, skipping merge")

    if config['RunSingleAnalysis']:
        for i, sample in enumerate(config["SampleIds"]):
            config_local = copy.deepcopy(config)
            config_local["AnalysisFolder"] = os.path.join(config["PathToAnalysis"],
                                                          config["SampleIds"][i])
            if config["Objects"]:
                config_local["Object"] = config["Objects"][i]
            else:
                config_local["SampleId"] = config["SampleIds"][i]
                config_local["PathToCount"] = config["PathsToCounts"][i]
            config_local["AnalysisType"] = "single"
            config_local["Organism"] = config["Organism"][i]
            if not config_local["FilterUMI"] or config_local['WholeUMI']:
                analysis_run(config_local)
                config_local['WholeUMI'] = False
            if config_local["FilterUMI"]:
                config_local["AnalysisFolder"] = os.path.join(config_local["AnalysisFolder"], "filtered")
                analysis_run(config_local)
