import argparse
import os
import subprocess


def generate_dirs(directory: str) -> None:
    if not directory.endswith('/'):
        directory = f'{directory}/'
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".yaml"):
            subprocess.call(['python', 'run.py', '--config', f'{directory}{filename}'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate directories with target pipeline")
    parser.add_argument('--yaml_path', type=str, required=True,
                        help='Path to the yaml configuration folder')
    args = parser.parse_args()
    generate_dirs(args.yaml_path)
