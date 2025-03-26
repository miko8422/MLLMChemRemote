#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2025/3/21 2:40
# @Author  : miko
# @File    : JSME_Gen.py
# @Software: PyCharm
# @Site:

# TODO: This is a 3D version.
"""import json
import urllib.request as urlreq
from dash import Dash, html, Input, Output, callback
import dash_bio as dashbio

app = Dash()


model_data = urlreq.urlopen(
    'https://git.io/mol2d_buckminsterfullerene.json'
).read().decode('utf-8')

model_data = json.loads(model_data)

app.layout = html.Div([
    dashbio.Molecule2dViewer(
        id='dashbio-default-molecule2d',
        modelData=model_data
    ),
    html.Hr(),
    html.Div(id='default-molecule2d-output')
])

@callback(
    Output('default-molecule2d-output', 'children'),
    Input('dashbio-default-molecule2d', 'selectedAtomIds')
)
def update_selected_atoms(ids):
    if ids is None or len(ids) == 0:
        return "No atom has been selected. Select atoms by clicking on them."
    return "Selected atom IDs: {}.".format(', '.join([str(i) for i in ids]))

if __name__ == '__main__':
    app.run(debug=True)"""


# TODO: 2D version
from dash import Dash, html
import dash_bio as dashbio

app = Dash()

app.layout = html.Div([
    dashbio.Jsme(
        options= "useOpenChemLib query exportSVG",
        smiles='O=C(Nc1cccc(Cl)c1)c3cncc4nnc(c2ccc(OC(F)F)cc2)n34',
    ),
])

if __name__ == '__main__':
    app.run(debug=True, port=7800)

