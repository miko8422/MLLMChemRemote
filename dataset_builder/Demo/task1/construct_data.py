#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2025/3/24 3:52
# @Author  : miko
# @File    : construct_data.py
# @Software: PyCharm
# @Site:

import json
import random

from dataset_builder.Demo.demo_utils.rdkit_smiles_recog import gen_smile_img


dataset_path = f"../../../dataset/train/train_sample.json"
img_path = './data_demo_task1/'

dataset = []
with open(f'./task1_demo.json', 'w', encoding='utf-8') as f_save:
    with open(dataset_path, 'r', encoding='utf-8') as f_read:
        smiles = json.loads(f_read.read())
        smiles = random.sample(smiles, 2)

        for line in smiles:
            try:
                gen_smile_img(line, file_path=img_path)
                dataset.append(
                    {
                        "messages": [
                            {
                                "role": "user",
                                "content": "You are a chemical assistant that outputs SMILES notation for molecular graphs. When shown a molecular graph image, respond ONLY with the correct SMILES formula string - no explanations or additional text. <image> What is the SMILES formula of this molecule?"
                            },
                            {
                                "role": "assistant",
                                "content": f"{line['smile']}"
                            }
                        ],
                        "images": [
                            f"{img_path}/{line['id']}.png"
                        ]
                    }
                )
            except Exception as e:
                print(f"Error Smile {deep_smiles}: {e}")

        if len(dataset) != len(smiles):
            pass
        else:
            f_save.write(json.dumps(dataset))


with open(f'./task1_demo.json', 'r', encoding='utf-8') as f:
    dataset = json.loads(f.read())
    print(len(dataset))
    print(type(dataset))
    print(dataset[0])