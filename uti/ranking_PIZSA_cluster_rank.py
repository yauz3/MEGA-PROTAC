# -*- coding: utf-8 -*-
import re
import os
import glob
import time
import argparse
import textwrap


def extract_z_score(filename):
    # Dosyayı aç ve satırları oku
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Z-score değerini bul
    z_score_line = [line for line in lines if 'Z-score' in line]
    if z_score_line:
        z_score_value = float(z_score_line[0].split(':')[1].strip())
        print("Z-score:", z_score_value)
    else:
        print("Z-score bulunamadı.")
    return z_score_value

def PIZSA_score(protein, input_path):
    # Geçerli çalışma dizinini al
    current_path = os.getcwd()

    # Üst dizini al
    parent_path = os.path.dirname(current_path)
    print("parent_path")
    print(parent_path)
    os.chdir("{}/bin/PIZSA".format(parent_path))

    print("python2 run_PIZSA.py {}/{} -o {}/outputs/Megadock_OUT/pizsa/{}".format(input_path,
            protein, parent_path, protein))
    os.system("python2 run_PIZSA.py {}/{} -o {}/outputs/Megadock_OUT/pizsa/{}".format(input_path,
            protein, parent_path, protein))

    os.chdir("{}/outputs/Megadock_OUT/pizsa".format(parent_path))
    filename=str(protein).replace("pdb","")
    output_files = glob.glob("*{}*scores.txt*".format(filename))[0]
    # output_files=glob.glob("*{}._score*".format(pro.split(".")[0]))[0]
    #output_files = "{}._scores.txt{}._scores.txt".format(file, file)
    print("output_files",output_files)
    z_score = extract_z_score(output_files)
    return z_score


def parseArguments():
    parser = argparse.ArgumentParser(prog='python2 ranking_PIZSA.py',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''\
                                     ------------------------------------------------------------------------
                                     >> The script is to make assessment of protein 3D structure.
                                     ========================================================================
                                     '''),
                                     epilog=textwrap.dedent('''\
                                     ========================================================================
                                     >>Please read instructions before start.
                                     ENJOY!!!
                                     ------------------------------------------------------------------------'''))
    # girdiğinin üzeirndekini değiştirmek için input tanımladık. terminalde -h olara çalıştır görürsün inputu
    input = parser.add_argument_group('------------------------------------------------------------------------'
                                      '\n>input')
    input.add_argument('-p', '--protein', help=textwrap.dedent("""List of proteins"""))
    input.add_argument('-i', '--input_path', help=textwrap.dedent("""protein path"""))

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parseArguments()
    print(PIZSA_score(args.protein, args.input_path))
