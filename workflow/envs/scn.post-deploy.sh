#!/usr/bin/env bash

set -euo pipefail

## lines below require GITHUB_PAT in ~/.Renviron
R -e "devtools::install_github(repo='ctlab/SCNPrep', ref='bef352ef6106ebe86bbb8f32dc23fffd6398e7a2', upgrade='never')"