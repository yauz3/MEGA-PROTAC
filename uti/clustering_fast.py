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


def cluster(input_path,output_name="cluster_5_3_model",min_number=3,threshold=0.5):
    os.chdir(input_path)
    print("first step")
    os.system(f"for pdb in $( cat pdb.list ); do python2 {FCC_path}/pdb_chainxseg.py $pdb > temp; mv temp $pdb; done")
    print("second step")
    os.system(f'python2 {FCC_path}/make_contacts.py -f pdb.list -n 8 -e {fcc}/src/contact_fcc')
    print("third step")
    os.system(r"sed -e 's/pdb/contacts/' pdb.list | sed -e '/^$/d' > pdb.contacts")
    print("fourth step")
    os.system(f'python2 {FCC_path}/calc_fcc_matrix.py -f pdb.contacts -o fcc_matrix.out')
    print("five step")
    os.system(f'python2 {FCC_path}/cluster_fcc.py fcc_matrix.out {threshold} -o clusters_{threshold}.out -c {min_number}')
    print("six step")
    os.system(
        f'python2 {FCC_path}/ppretty_clusters.py clusters_{threshold}.out pdb.list > {output_name}')
    print("seven ")
    return 'cluster_5_3_model'
