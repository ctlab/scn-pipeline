import os
import logging
import subprocess


def split_path_to_parts(path: str) -> list[str]:
    path = os.path.normpath(path)
    parts = path.split(os.sep)
    return parts


def join_path_parts(parts: list[str]) -> str:
    return os.path.join(*parts)


def box_upload_folder(folder_to_copy: str,
                      box_folder: str,
                      box_receipt: str):
    ptlf = os.path.normpath(folder_to_copy)
    ptl = os.path.normpath(box_folder)
    ptl_parts = split_path_to_parts(ptl)
    written_files = []

    for root, dirs, files in os.walk(ptlf):
        relative_path = root.replace(ptlf, "", 1).strip(os.sep)
        box_path = os.path.join(join_path_parts(ptl_parts), relative_path).strip(os.sep)
        for file in files:
            local_file_path = os.path.join(root, file)
            logging.warning(f"Uploading {local_file_path} to folder {box_path}")

            command = ["rclone", "copy", local_file_path, f"remote:{box_path}"]
            logging.warning(command)
            env = os.environ.copy()
            output = subprocess.check_output(" ".join(command), shell=True, env=env)
            written_files.append(local_file_path)

    with open(box_receipt, "w") as receipt:
        for wf in written_files:
            receipt.write(wf + "\n")


if __name__ == "__main__":
    folder_to_copy = snakemake.params["to_upload"]
    box_folder = snakemake.params["box_path"]
    box_receipt = snakemake.output["box_receipt"]
    box_upload_folder(folder_to_copy, box_folder, box_receipt)