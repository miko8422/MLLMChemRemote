#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2025/3/24 3:38
# @Author  : miko
# @File    : rdkit_smiles_recog.py
# @Software: PyCharm
# @Site:

import os
from rdkit import Chem
from rdkit.Chem import Draw
from deepsmiles import Converter
import matplotlib.pyplot as plt
import numpy as np

def gen_smile_img(smiles_info, file_path, size=(400, 300), show=False):
    """
        Generate one Smile Image
    :param smiles:
    :param file_path:
    :param size:
    :return:
    """
    file_path = file_path.rstrip('/')
    if os.path.isfile(file_path):
        raise ValueError(f"The path '{file_path}' is a file. Expected a directory.")
    elif not os.path.exists(file_path):
        os.makedirs(file_path)
        print(f"Created directory: {file_path}")

    smiles = smiles_info['smile']
    id = smiles_info['id']

    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol, size=size)
    plt.figure(figsize=(size[0] / 100, size[1] / 100))
    plt.imshow(np.array(img))
    plt.axis('off')
    plt.tight_layout()

    plt.savefig(f'{file_path}/{id}.png')
    if show:
        plt.show()
    else:
        plt.close()


if __name__ == '__main__':
    idx = 1
    smiles = "C1=CC=C(C=C1)N2C=C(C3=CC=CC=C32)/C=C(\\CN)/C(=O)O"

    smiles_infos = [{'id': idx, 'smile': smiles}]
    size = (400, 300)
    path = '../task1/data_demo_task1/'

    for line in smiles_infos:
        gen_smile_img(line, file_path=path)
