#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2025/3/24 23:29
# @Author  : miko
# @File    : read_demo.py
# @Software: PyCharm
# @Site:

# TODO: Here is just for demostration of the dataset.

import random
import json
from tabulate import tabulate

with open(f'task6/task6_demo.json', 'r', encoding='utf-8') as f:
    dataset = json.loads(f.read())

print(json.dumps(dataset, indent=2))