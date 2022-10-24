import sys
import os
import logging
import dropbox
from dropbox.files import WriteMode, UploadSessionCursor, CommitInfo
from dropbox.exceptions import ApiError, AuthError
from pathlib import Path


logging.basicConfig(level=logging.INFO)

DROPBOX_TOKEN = os.environ["DROPBOX_KEY"]
COMMON_PREFIX = "/ArtyomovLab Team Folder/scn-preloaded-dev-tmp/datasets/dump_2022/"
CHUNK_SIZE = 10 * 1024 * 1024


def upload_file_single_chunk(dbx: dropbox.Dropbox, local_file: str, file_location: str):
    file_stats = Path(local_file).stat()
    mb_size = file_stats.st_size / (1024 * 1024)
    with open(local_file, 'rb') as lf:
        logging.info(f"Uploading {local_file}, size - {mb_size:.2f}MB using file_upload")
        try:
            result = dbx.files_upload(lf.read(), file_location, mode=WriteMode('overwrite'))
            logging.info(f"Successful upload - {result.name}")
            return result
        except ApiError as err:
            if (err.error.is_path() and
                    err.error.get_path().reason.is_insufficient_space()):
                sys.exit("ERROR: Cannot back up; insufficient space.")
            elif err.user_message_text:
                print(err.user_message_text)
                sys.exit()
            else:
                print(err)
                sys.exit()


def upload_file_session(dbx: dropbox.Dropbox, local_file: str, file_location: str):
    file_stats = Path(local_file).stat()
    mb_size = file_stats.st_size / (1024 * 1024)
    with open(local_file, 'rb') as lf:
        logging.info(f"Uploading {local_file}, size - {mb_size:.2f}MB using upload_session with 10MB chunks")
        try:
            data = lf.read(CHUNK_SIZE)
            session = dbx.files_upload_session_start(data)
            session_id = session.session_id
            offset = len(data)
            cursor = UploadSessionCursor(session_id, offset)
            logging.info(f"Session {session_id} started, currently uploaded - {cursor.offset / (1024 * 1024):.2f}MB")

            data = lf.read(CHUNK_SIZE)
            while data:
                dbx.files_upload_session_append_v2(data, cursor)
                cursor.offset += len(data)
                logging.info(
                    f"Session {session_id} continued, currently uploaded - {cursor.offset / (1024 * 1024):.2f}MB")
                data = lf.read(CHUNK_SIZE)

            commit = CommitInfo(
                path=file_location,
                mode=WriteMode('overwrite')
            )
            result = dbx.files_upload_session_finish(bytes(), cursor, commit)
            logging.info(f"Session {session_id} finished - {result.name} uploaded")
            return result

        except ApiError as err:
            if (err.error.is_path() and
                    err.error.get_path().reason.is_insufficient_space()):
                sys.exit("ERROR: Cannot back up; insufficient space.")
            elif err.user_message_text:
                print(err.user_message_text)
                sys.exit()
            else:
                print(err)
                sys.exit()


def upload_file(dbx: dropbox.Dropbox, local_file: str, file_location: str):
    file_stats = Path(local_file).stat()
    # file_upload has a limit of 150M otherwise start an upload session
    if file_stats.st_size <= CHUNK_SIZE:
        result = upload_file_single_chunk(dbx, local_file, file_location)
    else:
        result = upload_file_session(dbx, local_file, file_location)
    return result


if __name__ == '__main__':
    if len(DROPBOX_TOKEN) == 0:
        sys.exit("ERROR: Looks like you didn't add your access token.")

    # Create an instance of a Dropbox class, which can make requests to the API.
    print("Creating a Dropbox object...")
    with dropbox.Dropbox(DROPBOX_TOKEN) as dbx:
        # Check that the access token is valid
        try:
            dbx.users_get_current_account()
        except AuthError:
            sys.exit("ERROR: Invalid access token; try re-generating an "
                     "access token from the app console on the web.")

    results_directory = os.path.dirname(snakemake.input["descriptor"])
    path_prefix = snakemake.params["path_prefix"]
    receipt = snakemake.output["dropbox_receipt"]


    with open(receipt, "w") as receipt:
        for root, subdirs, files in os.walk(results_directory):
            for file in files:
                path = os.path.join(root, file)
                relative_path = path.replace(results_directory, "").strip("/")
                resulting_path = os.path.join(COMMON_PREFIX, path_prefix, relative_path)
                result = upload_file(dbx, path, resulting_path)
                receipt.write(result.name + "\n")


