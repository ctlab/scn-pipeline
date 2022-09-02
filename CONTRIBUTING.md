# Contributing

When contributing to this repository, please first discuss the change you wish to make via issue,
email, or any other method with the owners of this repository before making a change. 

Please note we have a code of conduct, please follow it in all your interactions with the project.

## Main branch

While this repository is used to generate files for scNavigator, the scope of the main branch
is smaller: get a processed Seurat object from a ID in a public database (like GEO). We believe that
smaller scope of the main branch will allow more people to use this repository for their own purposes.
Please, keep it in main when contribute to main branch.

## Before creating a pull request

1. Before contributing to this repository, please create an issue describing things 
that you intend to solve with your changes to this repository. We might be already working on these,
or we might indeed welcome your input on solving the issue.
2. We consider this repository to be a "Snakemake pipeline" so please refer to 
[snakemake documentation](https://snakemake.readthedocs.io/en/stable/) for documentation and best practices.
3. When creating a new `rule` in this pipeline, please create **minimal** [conda](https://docs.conda.io/en/latest/) 
environment required to run this rule. If the rule uses a tool already used in the pipeline, 
please reuse the existing environment.
4. Please make sure your changes don't break any tests. If you think that the test in incorrect, please,
report this in a separate issue.

## Pull request process

1. Please use [conventional commits](https://www.conventionalcommits.org/en/v1.0.0/) (this is required for releases) and provide context to your
changes and reference the issue this commit is solving.
2. We use [release-please](https://github.com/googleapis/release-please) bot to make releases and keep versioning
scheme of [SemVer](http://semver.org/).
3. You may merge the Pull Request in once you have the sign-off of other developers (@konsolerr), or if you 
   do not have permission to do that, you may request the developers (@konsolerr) to merge it for you.