import pandas as pd
import os


def txt_to_csv(txt_file_path):
    txtin = open(txt_file_path, "rt")
    out_path = txt_file_path.replace(".txt", ".csv")
    csvout = open(out_path, "wt")

    for line in txtin:
        csvout.write(' '.join(line.split()))
        csvout.write(' \n')

    txtin.close()
    csvout.close()

    readFile = pd.read_csv(txt_file_path,
                           delim_whitespace=True,  on_bad_lines='skip', encoding='unicode_escape')

    readFile.to_csv(out_path, index=None)

    """This will be used to append to a list obj in R, so that it can iterate through that
    path list and perform another process"""

    return out_path


def txt_to_csv_files_in_root(root_path):
    sumstat_out_paths = []
    # The return serves as a list the next process will need, which is in R
    for root, dirs, files in os.walk(root_path):
        for name in files:
            if name.endswith(".txt") and not name.startswith("._"):
                file_path = os.path.join(root, name)
                sumstat_out_path = txt_to_csv(file_path)
                sumstat_out_paths.append(sumstat_out_path)

    return sumstat_out_paths
