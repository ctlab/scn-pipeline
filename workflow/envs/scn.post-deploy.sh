#!/usr/bin/env bash

set -euo pipefail

## lines below require GITHUB_PAT in ~/.Renviron
R -e "devtools::install_github(repo='ctlab/SCNPrep', ref='427b87b52e7fbc07c056aeb4e8729f361b075745', upgrade='never')"