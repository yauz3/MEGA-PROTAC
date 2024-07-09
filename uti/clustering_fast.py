import glob
import os,re
import shutil
import re
import subprocess


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)) # ..
FCC_path = os.path.join(PROJECT_ROOT, "bin/fcc/scripts")
voronota_path = os.path.join(PROJECT_ROOT, "bin/voronota")
fcc= os.path.join(PROJECT_ROOT, "bin/fcc")

def prepare_input_file(input_path,file_list):
    os.chdir(input_path)
    with open("pdb.list", "w") as pdb_list_file:
        for file_name in file_list:
            pdb_list_file.write(f"{file_name}.pdb" + "\n")

    return "pdb.list"


def cluster(input_path, output_name="cluster_5_3_model", min_number=3, threshold=0.5):
    os.chdir(input_path)

    print("First step")
    first_command = f"for pdb in $(cat pdb.list); do python2 {FCC_path}/pdb_chainxseg.py $pdb > temp; mv temp $pdb; done"
    subprocess.run(["conda", "run", "-n", "pizsa", "bash", "-c", first_command], capture_output=True, text=True)

    print("Second step")
    second_command = [
        "conda", "run", "-n", "pizsa", "python2", f"{FCC_path}/make_contacts.py", "-f", "pdb.list", "-n", "8", "-e",
        f"{fcc}/src/contact_fcc"
    ]
    subprocess.run(second_command, capture_output=True, text=True)

    print("Third step")
    third_command = "sed -e 's/pdb/contacts/' pdb.list | sed -e '/^$/d' > pdb.contacts"
    subprocess.run(["conda", "run", "-n", "pizsa", "bash", "-c", third_command], capture_output=True, text=True)

    print("Fourth step")
    fourth_command = [
        "conda", "run", "-n", "pizsa", "python2", f"{FCC_path}/calc_fcc_matrix.py", "-f", "pdb.contacts", "-o",
        "fcc_matrix.out"
    ]
    subprocess.run(fourth_command, capture_output=True, text=True)

    print("Fifth step")
    fifth_command = [
        "conda", "run", "-n", "pizsa", "python2", f"{FCC_path}/cluster_fcc.py", "fcc_matrix.out", str(threshold), "-o",
        f"clusters_{threshold}.out", "-c", str(min_number)
    ]
    subprocess.run(fifth_command, capture_output=True, text=True)

    print("Sixth step")
    
    sixth_command = f"conda run -n pizsa python2 {FCC_path}/ppretty_clusters.py clusters_{threshold}.out pdb.list > {output_name}"
    subprocess.run(["bash", "-c", sixth_command], capture_output=True, text=True)

    print("Clustering is finished!")
    return 'cluster_5_3_model'
