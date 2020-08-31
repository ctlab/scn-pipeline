import yaml

from analysis_run_merge import analysis


def pipeline(config_file_name: str) -> None:
    with open(config_file_name) as config, \
            open("defaults.yaml") as defaults:
        run_dict = yaml.load(defaults, Loader=yaml.FullLoader)
        config_dict = yaml.load(config, Loader=yaml.FullLoader)
        run_dict.update(config_dict)
    if type(run_dict['Objects']) is str:
        assert type(run_dict['SampleIds']) == type(
            run_dict['Objects']), "Probably you specify one SampleIds and several Objects or vice versa"
        run_dict['Objects'] = run_dict['Objects'].split()
        run_dict['SampleIds'] = run_dict['SampleIds'].split()
    assert run_dict['RunName'], "You must specify RunName"
    assert run_dict['SampleIds'], "Number of SampleIds can not be zero"
    assert len(run_dict['SampleIds']) == len(run_dict['Objects']), "SampleIds number must correspond to Objects number"
    assert (run_dict['Objects'] and not run_dict['PathsToCounts']) or (
            not run_dict['Objects'] and run_dict['PathsToCounts']), "Specify Objects or PathsToCounts"
    assert run_dict['PathToAnalysis'], "PathToAnalysis can not be empty"
    analysis(run_dict)

