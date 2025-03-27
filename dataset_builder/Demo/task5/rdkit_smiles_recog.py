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
import os
from PIL import Image, ImageDraw, ImageFont

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

def gen_multiple_smile_img(smiles_list, file_path, id, size=(400, 300), show=False, add_arrow=True, add_question=True):
    """
    Generate a single image containing multiple molecules with '+' symbols between them,
    followed by an arrow and optionally a question mark

    :param smiles_list: List of SMILES strings
    :param file_path: Path to save the image
    :param id: ID for the file name
    :param size: Base size for each molecule image
    :param show: Whether to show the image
    :param add_arrow: Whether to add an arrow at the end
    :param add_question: Whether to add a question mark after the arrow
    :return:
    """
    file_path = file_path.rstrip('/')
    if os.path.isfile(file_path):
        raise ValueError(f"The path '{file_path}' is a file. Expected a directory.")
    elif not os.path.exists(file_path):
        os.makedirs(file_path)
        print(f"Created directory: {file_path}")

    # Generate molecules and images
    mols = [Chem.MolFromSmiles(smile) for smile in smiles_list]
    mol_images = [Draw.MolToImage(mol, size=size) for mol in mols]

    # Define symbol widths
    plus_width = 40  # Space for '+' symbol
    arrow_width = 80  # Space for arrow symbol
    question_width = 40  # Space for question mark

    # Calculate total width
    symbols_width = (len(mol_images) - 1) * plus_width
    if add_arrow:
        symbols_width += arrow_width
    if add_question:
        symbols_width += question_width

    total_width = sum([img.width for img in mol_images]) + symbols_width
    max_height = max([img.height for img in mol_images])

    # Create a new combined image
    combined_img = Image.new('RGB', (total_width, max_height), (255, 255, 255))
    draw = ImageDraw.Draw(combined_img)

    # Try to load a font, use default if not available
    try:
        font = ImageFont.truetype("arial.ttf", 40)
    except IOError:
        # If arial isn't available, try a system font or load default
        try:
            font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 40)
        except IOError:
            font = ImageFont.load_default()

    # Paste all images and add symbols
    x_offset = 0
    for i, img in enumerate(mol_images):
        # Paste the molecule image
        combined_img.paste(img, (x_offset, (max_height - img.height) // 2))
        x_offset += img.width

        # Add '+' symbol if not the last molecule
        if i < len(mol_images) - 1:
            draw.text((x_offset + 10, max_height // 2 - 20), "+", fill=(0, 0, 0), font=font)
            x_offset += plus_width

    # Add arrow symbol after all molecules
    if add_arrow:
        draw.text((x_offset + 10, max_height // 2 - 20), "→", fill=(0, 0, 0), font=font)
        x_offset += arrow_width

    # Add question mark after arrow
    if add_question:
        draw.text((x_offset, max_height // 2 - 20), "?", fill=(0, 0, 0), font=font)

    # Save the combined image
    combined_img.save(f'{file_path}/{id}.png')

    if show:
        plt.figure(figsize=(total_width / 100, max_height / 100))
        plt.imshow(np.array(combined_img))
        plt.axis('off')
        plt.tight_layout()
        plt.show()

    return combined_img


def gen_reaction_img(reactants_smiles, product_smiles, file_path, id, size=(400, 300), show=False):
    """
    Generate a reaction image with reactants connected by + symbols and an arrow to the product

    :param reactants_smiles: List of SMILES strings for reactants
    :param product_smiles: SMILES string for product
    :param file_path: Path to save the image
    :param id: ID for the file name
    :param size: Base size for each molecule image
    :param show: Whether to show the image
    :return: The combined image
    """
    file_path = file_path.rstrip('/')
    if os.path.isfile(file_path):
        raise ValueError(f"The path '{file_path}' is a file. Expected a directory.")
    elif not os.path.exists(file_path):
        os.makedirs(file_path)
        print(f"Created directory: {file_path}")

    # Generate molecules and images for reactants
    reactant_mols = []
    for smile in reactants_smiles:
        mol = Chem.MolFromSmiles(smile)
        if mol is not None:
            reactant_mols.append(mol)
        else:
            print(f"Warning: Could not create molecule from SMILES: {smile}")

    if not reactant_mols:
        print(f"Error: No valid reactant molecules could be created for ID {id}")
        return None

    # Generate molecule for product
    product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

    # Create images
    reactant_images = [Draw.MolToImage(mol, size=size) for mol in reactant_mols]
    product_image = Draw.MolToImage(product_mol, size=size) if product_mol else None

    # Define symbol widths
    plus_width = 40  # Space for '+' symbol
    arrow_width = 80  # Space for arrow symbol

    # Calculate total width
    symbols_width = (len(reactant_images) - 1) * plus_width + arrow_width

    total_width = sum([img.width for img in reactant_images])
    if product_image:
        total_width += product_image.width
    total_width += symbols_width

    # Determine max height
    all_images = reactant_images + ([product_image] if product_image else [])
    max_height = max([img.height for img in all_images]) if all_images else size[1]

    # Create a new combined image
    combined_img = Image.new('RGB', (total_width, max_height), (255, 255, 255))
    draw = ImageDraw.Draw(combined_img)

    # Try to load a font, use default if not available
    try:
        font = ImageFont.truetype("arial.ttf", 40)
    except IOError:
        # If arial isn't available, try a system font or load default
        try:
            font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 40)
        except IOError:
            font = ImageFont.load_default()

    # Paste reactant images and add '+' symbols
    x_offset = 0
    for i, img in enumerate(reactant_images):
        # Paste the molecule image
        combined_img.paste(img, (x_offset, (max_height - img.height) // 2))
        x_offset += img.width

        # Add '+' symbol if not the last reactant
        if i < len(reactant_images) - 1:
            draw.text((x_offset + 10, max_height // 2 - 20), "+", fill=(0, 0, 0), font=font)
            x_offset += plus_width

    # Add arrow symbol after all reactants
    draw.text((x_offset + 10, max_height // 2 - 20), "→", fill=(0, 0, 0), font=font)
    x_offset += arrow_width

    # Add product image
    if product_image:
        combined_img.paste(product_image, (x_offset, (max_height - product_image.height) // 2))
    else:
        # If no product, add a question mark
        draw.text((x_offset, max_height // 2 - 20), "?", fill=(0, 0, 0), font=font)

    # Save the combined image
    combined_img.save(f'{file_path}/{id}.png')

    if show:
        plt.figure(figsize=(total_width / 100, max_height / 100))
        plt.imshow(np.array(combined_img))
        plt.axis('off')
        plt.tight_layout()
        plt.show()

    return combined_img

if __name__ == '__main__':
    idx = 1
    # smiles = "C1=CC=C(C=C1)N2C=C(C3=CC=CC=C32)/C=C(\\CN)/C(=O)O"

    smiles_list = [
        "C1=CC=C(C=C1)N2C=C(C3=CC=CC=C32)/C=C(\\CN)/C(=O)O",
        "CC(=O)O",  # Add your additional reactants here
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    ]

    smiles_infos = [{'id': idx, 'smile': smiles} for idx, smiles in enumerate(smiles_list)]
    size = (400, 300)
    path = './data_demo_task3/'

    gen_reaction_img(smiles_list, file_path=path, id=idx, size=size)


    # for line in smiles_infos:
    #     gen_smile_img(line, file_path=path)
