import os

import urllib.request



# Folder to save three data files

data_dir = os.path.expanduser("~/capstone_epiclock_project/data_files")

if not os.path.exists(data_dir):

    os.makedirs(data_dir)



# FTP links of files from GEO

files = {

    "metadata": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87571/suppl/GSE87571%5Fadditional%5Fsample%5Fchararcteristics.xlsx",

    "beta_matrix1": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87571/suppl/GSE87571%5Fmatrix1of2.txt.gz",

    "beta_matrix2": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87571/suppl/GSE87571%5Fmatrix2of2.txt.gz"

}



# Download each file

for name, url in files.items():

    output_path = os.path.join(data_dir, url.split("/")[-1])

    print(f"Downloading {name}...")

    urllib.request.urlretrieve(url, output_path)

    print(f"Saved to {output_path}")



print("All three files downloaded!")

