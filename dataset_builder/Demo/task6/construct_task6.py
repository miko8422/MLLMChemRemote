#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2025/3/24 4:03
# @Author  : miko
# @File    : construct_task2.py
# @Software: PyCharm
# @Site:

import json
import random
from datasets import load_dataset
import selfies as sf

from rdkit_smiles_recog import gen_smile_img


def selfies2smiles(selfies):
    return sf.decoder(selfies)


dataset = load_dataset("zjunlp/Mol-Instructions", "Molecule-oriented Instructions", trust_remote_code=True)
mol_describe = dataset['property_prediction']
mol_describe = mol_describe.shuffle(seed=42)
print(mol_describe[0])
# print(mol_describe[0]['input'].split('.'))
# raise Exception

dataset_path = f"../../../dataset/train/train_sample.json"
img_path = './data_demo_task6/'

dataset = []
with open(f'./task6_demo.json', 'w', encoding='utf-8') as f_save:
    two_examples = mol_describe.select(range(2))
    smiles = [{'id': idx, 'smile': selfies2smiles(line['input']), 'input': line['input'], 'output': line['output']} for
              idx, line in enumerate(two_examples)]
    for line in smiles:
        try:
            gen_smile_img(line, file_path=img_path)
            dataset.append(
                {
                    "messages": [
                        {
                            "role": "user",
                            "content": "You are a chemical assistant that analyzes molecular structures and provides precise calculations for HOMO-LUMO gap energy. <image> I need to know the HOMO-LUMO gap energy of this molecule, could you please provide it?"
                        },
                        {
                            "role": "assistant",
                            "content": f"{line['output']}"
                        }
                    ],
                    "images": [
                        f"{img_path}/{line['id']}.png"
                    ]
                }
            )
        except Exception as e:
            print(f"Error Smile {line}: {e}")

    if len(dataset) != len(smiles):
        pass
    else:
        f_save.write(json.dumps(dataset))


with open(f'./task6_demo.json', 'r', encoding='utf-8') as f:
    dataset = json.loads(f.read())
    print(len(dataset))
    print(type(dataset))
    print(dataset[0])


# # Test Demo
# if __name__ == "__main__":
#     dataset = load_dataset("zjunlp/Mol-Instructions", "Molecule-oriented Instructions", trust_remote_code=True)
#     mol_describe = dataset['molecular_description_generation']
#     sample_data = mol_describe.shuffle(seed=42)[0]
#     print(sample_data['input'])
#     print(sf.decoder(sample_data['input']))


