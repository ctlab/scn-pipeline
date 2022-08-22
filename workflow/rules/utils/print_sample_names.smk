from workflow.scripts.DependencyDispatcher import DependencyDispatcher

rule print_sample_names:
  input: config['samples']
  run:
    datasets = DependencyDispatcher(config).get_all_datasets()
    for dataset in datasets:
      wildcards = DependencyDispatcher.get_wildcards_placeholder()
      wildcards.dataset = dataset
      print("######")
      print(f"Showing sample info for dataset {dataset}")
      print(DependencyDispatcher(config).get_sample_names(wildcards))