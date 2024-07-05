# -*- coding: utf-8 -*-
# 28/06/2024
# Author: Sadettin Y. Ugurlu & David McDonald


import re


# Setting the seed for the random number generator to ensure reproducibility
# The seed value (42) is arbitrary but consistent, so the sequence of random numbers
# generated will be the same each time this code is run
random.seed(42)

# Get the current working directory
current_path = os.getcwd()


def group_values(protein,selected_protein):
    re_cluster = re.compile('Cluster (.*?) -> (.*?)\n', re.S)
    re_number = re.compile(f'(.*?)_{protein}_target_input_pro_comp_(\d+)')
    with open(f"{current_path}/outputs/Megadock_OUT/{protein}_rotate/cluster_5_3_model_filtered", 'r') as cluster_file:
        cluster_lines = cluster_file.readlines()

    for cluster_line in cluster_lines:
        re_line = re.search(re_cluster, cluster_line)
        if re_line:
            cluster_id, model = re_line.group(1), re_line.group(2)
            model_list = model.split()
            if selected_protein in model_list:
                print(model_list)

group_values("7q2j-CD",
             "7q2j-CD_target_input_0_3895_center_x_0-0_y_-3-0_z_4-5_center_x_0-0_y_0-0_z_-10-0_center"
             )