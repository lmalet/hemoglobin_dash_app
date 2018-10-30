# -*- coding: utf-8 -*-
import dash
from dash.dependencies import Input, State, Output, Event
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import pandas as pd
import os
from Bio import SeqIO
from Bio.SeqUtils import ProtParam as pp
from Bio.SeqUtils import seq3
from plotly import tools

app = dash.Dash()
server = app.server

##### 1. DATA #####

    ##### 1.1 STRUCTURE #####
df1 = pd.read_csv("data/1a3n.csv")
df2 = pd.read_csv("data/2gdm.csv")
df3 = pd.read_csv("data/2qsp.csv")
df4 = pd.read_csv("data/3hrw.csv")
df5 = pd.read_csv("data/4odc.csv")

    ##### 1.2 COMPOSITION #####
fasta_1a3n = str(os.path.join("data","1a3n.fasta"))
fasta_2gdm = str(os.path.join("data","2gdm.fasta"))
fasta_2qsp = str(os.path.join("data","2qsp.fasta"))
fasta_3hrw = str(os.path.join("data","3hrw.fasta"))
fasta_4odc = str(os.path.join("data","4odc.fasta"))

##### 2. OPTIONS #####
option_species = {"Homo sapiens (Human)": "v1a3n",
                  "Mus musculus (House mouse)":"v3hrw",
                  "Trematomus bernacchii (Emerald rockcod)":"v4odc",
                  "Bos taurus (Cattle)":"v2qsp",
                  "Lupinus luteus (Yellow lupin)":"v2gdm"
                  }

##### 3. FUNCTIONS IN CALLBACKS #####

    ##### 3.1 for COMPOSITION #####
def choose_fasta(value):
    """
    choose the right fasta file
    :param value:
    :return:
    """
    switcher = {
        "v1a3n": fasta_1a3n,
        "v2gdm": fasta_2gdm,
        "v2qsp": fasta_2qsp,
        "v3hrw": fasta_3hrw,
        "v4odc": fasta_4odc
    }
    return switcher.get(value)

def get_count_aa(value):
    traces = []
    file_path = choose_fasta(value)
    with open(file_path, "r") as file_fasta:
        for entry in SeqIO.parse(file_fasta, "fasta"):
            id_prot = entry.id.split("|")
            id_chain = id_prot[0].split(":")
            seq = str(entry.seq)
            X = pp.ProteinAnalysis(seq)
            count_aa = X.count_amino_acids()
            aa_list = []
            aa_count = []
            for key, value in count_aa.items():
                aa_name = seq3(key)
                aa_list.append(aa_name)
                aa_count.append(value)
            traces.append(go.Bar(
                x=aa_list,
                y=aa_count,
                name="Chain " + id_chain[1]
                ))
    return traces


def get_percentage_aa(value):
    traces = []
    file_path = choose_fasta(value)
    with open(file_path, "r") as file_fasta:
        for entry2 in SeqIO.parse(file_fasta, "fasta"):
            id_prot = entry2.id.split("|")
            id_chain = id_prot[0].split(":")
            seq = str(entry2.seq)
            X = pp.ProteinAnalysis(seq)
            percent_aa = X.get_amino_acids_percent()
            aa_list = []
            aa_percent = []
            for key, value in percent_aa.items():
                aa_name = seq3(key)
                aa_list.append(aa_name)
                aa_percent.append(value*100)
            traces.append(go.Bar(
                x=aa_list,
                y=aa_percent,
                name="Chain " + id_chain[1]
                ))
    return traces


def get_weight(value):
    traces = []
    file_path = choose_fasta(value)
    with open(file_path, "r") as file_fasta:
        for entry in SeqIO.parse(file_fasta, "fasta"):
            id_prot = entry.id.split("|")
            id_chain = id_prot[0].split(":")
            seq = str(entry.seq)
            X = pp.ProteinAnalysis(seq)
            weight = X.molecular_weight()
            chain_name = []
            chain_weight = []
            chain_name.append("Chain " + id_chain[1])
            chain_weight.append(weight)
            traces.append(go.Bar(
                x=chain_name,
                y=chain_weight,
                name="Chain " + id_chain[1]
                ))
    return traces


count_1a3n = get_count_aa("v1a3n")
count_2gdm = get_count_aa("v2gdm")
count_2qsp = get_count_aa("v2qsp")
count_3hrw = get_count_aa("v3hrw")
count_4odc = get_count_aa("v4odc")


def choose_count(value):
    switcher = {
        "v1a3n": count_1a3n,
        "v2gdm": count_2gdm,
        "v2qsp": count_2qsp,
        "v3hrw": count_3hrw,
        "v4odc": count_4odc
    }
    return switcher.get(value)


percentage_1a3n = get_percentage_aa("v1a3n")
percentage_2gdm = get_percentage_aa("v2gdm")
percentage_2qsp = get_percentage_aa("v2qsp")
percentage_3hrw = get_percentage_aa("v3hrw")
percentage_4odc = get_percentage_aa("v4odc")


def choose_percentage(value):
    switcher = {
        "v1a3n": percentage_1a3n,
        "v2gdm": percentage_2gdm,
        "v2qsp": percentage_2qsp,
        "v3hrw": percentage_3hrw,
        "v4odc": percentage_4odc
    }
    return switcher.get(value)


weight_1a3n = get_weight("v1a3n")
weight_2gdm = get_weight("v2gdm")
weight_2qsp = get_weight("v2qsp")
weight_3hrw = get_weight("v3hrw")
weight_4odc = get_weight("v4odc")


def choose_weight(value):
    switcher = {
        "v1a3n": weight_1a3n,
        "v2gdm": weight_2gdm,
        "v2qsp": weight_2qsp,
        "v3hrw": weight_3hrw,
        "v4odc": weight_4odc
    }
    return switcher.get(value)


    ##### 3.2 for STRUCTURE #####
def color_picker(atom):
    switcher = {
        "N": "rgba(0, 0, 255, 0.8)",
        "C": "rgba(0, 0, 0, 0.8)",
        "O": "rgba(255, 0, 0, 0.8)",
        "S": "rgba(255, 128, 0, 0.8)"
    }
    return switcher.get(atom)


def hem_atom_color_picker(atom):
    switcher = {
        "N": "rgba(0, 0, 255, 1.0)",
        "C": "rgba(0, 0, 0, 1.0)",
        "O": "rgba(255, 0, 0, 1.0)",
        "FE": "rgba(44, 160, 44, 1.0)"
    }
    return switcher.get(atom)


def aa_color_picker(residue):
    switcher = {
        "ALA": "rgb(31,119,180)",
        "ARG": "rgb(255,127,14)",
        "ASN": "rgb(44,160,44)",
        "ASP": "rgb(214,39,40)",
        "CYS": "rgb(148,103,189)",
        "GLN": "rgb(140,86,75)",
        "GLU": "rgb(227,119,194)",
        "GLY": "rgb(127,127,127)",
        "HIS": "rgb(188,189,34)",
        "ILE": "rgb(23,190,207)",
        "LEU": "rgba(31,119,180,0.8)",
        "LYS": "rgba(255,127,14,0.8)",
        "MET": "rgba(44,160,44,0.8)",
        "PHE": "rgba(214,39,40,0.8)",
        "PRO": "rgba(148,103,189,0.8)",
        "PYL": "rgba(140,86,75,0.8)",
        "SER": "rgba(227,119,194,0.8)",
        "SEC": "rgba(127,127,127,0.8)",
        "THR": "rgba(188,189,34,0.8)",
        "TRP": "rgba(23,190,207,0.8)",
        "TYR": "rgba(31,119,180,0.6)",
        "VAL": "rgba(255,127,14,0.6)",
        "ASX": "rgba(44,160,44,0.6)",
        "GLX": "rgba(214,39,40,0.6)",
        "HEM": "rgb(255,0,0)"
    }
    return switcher.get(residue)



app.layout = html.Div([
    html.Div([
        html.H2(
            'Structure Analysis of the Hemoglobin Protein in Different Species',
            id='title'
        ),
        # html.Img(
        #     src="https://s3-us-west-1.amazonaws.com/plotly-tutorials/logo/new-branding/dash-logo-by-plotly-stripe-inverted.png"
        # )
    ],
        className="banner"
    ),
    html.Div([
        html.H5(
            "Choose the species"
        ),
        dcc.Dropdown(
            id="species-name",
            options=[
                {"label": key, "value": value} for key,value in option_species.items()
            ],
            value="v1a3n"
        )
    ],
        style={"padding": "20px"}
    ),
        html.Div([
            html.H5(
                "Choose the structure",
            ),
            dcc.Dropdown(
                id="structure-type",
                options=[
                    {"label": "Amino Acids", "value": "AA"},
                    {"label": "Atoms", "value": "atoms"},
                    {"label": "Subunits", "value": "subunits"}
                ],
                value="AA",
            ),
            html.Div([
                dcc.Graph(id="3D-structure")
            ],
                # style={"width": "48%"}
            )
        ],
            style={
                "width": "48%",
                # "padding-left": "20px",
                # "padding-right": "20px",
                "margin-left": "20px",
                "margin-right": "20px",
                "float": "left",
                "display": "inline-block"
            }
        ),
        html.Div([
            html.H5(
                "Composition"
            ),
            dcc.RadioItems(
                id="composition-radio",
                options=[
                    {"label": "Amino acids count", "value": "count"},
                    {"label": "Amino acids percentage", "value": "percentage"},
                    {"label": "Molecular weight", "value": "weight"}
                ],
                value="count",
                labelStyle={"display": "inline-block"}
            ),
            html.Div(
                id="first-species-selected"
            ),
            html.P("Compare with:"),

            dcc.RadioItems(
                id="second-species",
                options=[
                    {"label": "Homo sapiens", "value": "v1a3n"},
                    {"label": "Mus musculus", "value": "v3hrw"},
                    {"label": "Trematomus bernacchii", "value": "v4odc"},
                    {"label": "Bos taurus", "value": "v2qsp"},
                    {"label": "Lupinus luteus", "value": "v2gdm"}
                ],
                value="v3hrw",
                labelStyle={'display': 'inline-block'}
            ),
            dcc.Graph(id="composition-graph")
        ],
            style={
                # "width": "45%",
                "padding-right": "20px",
                "padding-left": "20px",
                "float": "left",
                "display": "inline-block"
            }
        ),
        html.Div([
            html.H6("References:"),
            html.P(["Data come from the ",
                    html.A("Protein Data Bank", href="https://dash.plot.ly/dash-html-components")
                    ]),
            html.Ul([
                html.Li("1A3N: Deoxy Human Hemoglobin"),
                html.Li("2GDM: Leghemoglobin (OXY)"),
                html.Li("2QSP: Bovine Hemoglobin at pH 5.7"),
                html.Li("3HRW: Crystal structure of hemoglobin from mouse (Mus musculus)at 2.8"),
                html.Li("4ODC: Crystal structure of Trematomus bernacchii hemoglobin in a partially cyanided state"),
            ]

            )
        ],
            style={"clear":"both",
                   # "display":"block",
                   "padding-left":"20px"}
        )
])


external_css = [
    # "https://cdnjs.cloudflare.com/ajax/libs/normalize/7.0.0/normalize.min.css",  # Normalize the CSS
    # "https://fonts.googleapis.com/css?family=Open+Sans|Roboto"  # Fonts
    # "https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css",
    # "https://cdn.rawgit.com/TahiriNadia/styles/faf8c1c3/stylesheet.css",
    # "https://cdn.rawgit.com/TahiriNadia/styles/b1026938/custum-styles_phyloapp.css"
]

for css in external_css:
    app.css.append_css({"external_url": css})


@app.callback(
    Output("3D-structure", "figure"),
    [Input("species-name", "value"),
     Input("structure-type", "value")]
)
def update_graph(pdb_value, value_type):
    ##### AMINO ACIDS #####
    if pdb_value == "v1a3n" and value_type == "AA":
        traces = []
        for residue in df1[df1["type"] == "ATOM"]["residue_name"].unique():
            traces.append(go.Scatter3d(
                x=df1[df1["residue_name"] == residue]["x_coordinate"],
                y=df1[df1["residue_name"] == residue]["y_coordinate"],
                z=df1[df1["residue_name"] == residue]["z_coordinate"],
                name=residue,
                mode="markers",
                marker=dict(
                    color=aa_color_picker(residue),
                    line=dict(width=1)
                )
            ))

        traces.append(go.Scatter3d(
                x=df1[df1["residue_name"] == "HEM"]["x_coordinate"],
                y=df1[df1["residue_name"] == "HEM"]["y_coordinate"],
                z=df1[df1["residue_name"] == "HEM"]["z_coordinate"],
                name="Heme",
                mode="markers",
                marker=dict(
                    color=aa_color_picker("HEM"),
                    line=dict(width=1)
                )
        ))
        return {
            "data": traces,
            "layout": go.Layout(
                height=700,
            )
        }
    elif pdb_value == "v2gdm" and value_type == "AA":
        traces = []
        for residue in df2[df2["type"] == "ATOM"]["residue_name"].unique():
            traces.append(go.Scatter3d(
                x=df2[df2["residue_name"] == residue]["x_coordinate"],
                y=df2[df2["residue_name"] == residue]["y_coordinate"],
                z=df2[df2["residue_name"] == residue]["z_coordinate"],
                name=residue,
                mode="markers",
                marker=dict(
                    color=aa_color_picker(residue),
                    line=dict(width=1)
                )
            ))

        traces.append(go.Scatter3d(
                x=df2[df2["residue_name"] == "HEM"]["x_coordinate"],
                y=df2[df2["residue_name"] == "HEM"]["y_coordinate"],
                z=df2[df2["residue_name"] == "HEM"]["z_coordinate"],
                name="Heme",
                mode="markers",
                marker=dict(
                    color=aa_color_picker("HEM"),
                    line=dict(width=1)
                )
        ))
        return {
            "data": traces,
            "layout": go.Layout(
                height=700,
            )
        }
    elif pdb_value == "v2qsp" and value_type == "AA":
        traces = []
        for residue in df3[df3["type"] == "ATOM"]["residue_name"].unique():
            traces.append(go.Scatter3d(
                x=df3[df3["residue_name"] == residue]["x_coordinate"],
                y=df3[df3["residue_name"] == residue]["y_coordinate"],
                z=df3[df3["residue_name"] == residue]["z_coordinate"],
                name=residue,
                mode="markers",
                marker=dict(
                    color=aa_color_picker(residue),
                    line=dict(
                        width=1
                    )
                )
            ))

        traces.append(go.Scatter3d(
                x=df3[df3["residue_name"] == "HEM"]["x_coordinate"],
                y=df3[df3["residue_name"] == "HEM"]["y_coordinate"],
                z=df3[df3["residue_name"] == "HEM"]["z_coordinate"],
                name="Heme",
                mode="markers",
                marker=dict(
                    color=aa_color_picker("HEM"),
                    line=dict(width=1)
                )
        ))
        return {
            "data": traces,
            "layout": go.Layout(
                height=700,
            )
        }
    elif pdb_value == "v3hrw" and value_type == "AA":
        traces = []
        for residue in df4[df4["type"] == "ATOM"]["residue_name"].unique():
            traces.append(go.Scatter3d(
                x=df4[df4["residue_name"] == residue]["x_coordinate"],
                y=df4[df4["residue_name"] == residue]["y_coordinate"],
                z=df4[df4["residue_name"] == residue]["z_coordinate"],
                name=residue,
                mode="markers",
                marker=dict(
                    color=aa_color_picker(residue),
                    line=dict(
                        width=1
                    )
                )
            ))

        traces.append(go.Scatter3d(
                x=df4[df4["residue_name"] == "HEM"]["x_coordinate"],
                y=df4[df4["residue_name"] == "HEM"]["y_coordinate"],
                z=df4[df4["residue_name"] == "HEM"]["z_coordinate"],
                name="Heme",
                mode="markers",
                marker=dict(
                    color=aa_color_picker("HEM"),
                    line=dict(width=1)
                )
        ))
        return {
            "data": traces,
            "layout": go.Layout(
                height=700,
            )
        }
    elif pdb_value == "v4odc" and value_type == "AA":
        traces = []
        for residue in df5[df5["type"] == "ATOM"]["residue_name"].unique():
            traces.append(go.Scatter3d(
                x=df5[df5["residue_name"] == residue]["x_coordinate"],
                y=df5[df5["residue_name"] == residue]["y_coordinate"],
                z=df5[df5["residue_name"] == residue]["z_coordinate"],
                name=residue,
                mode="markers",
                marker=dict(
                    color=aa_color_picker(residue),
                    line=dict(
                        width=1
                    )
                )
            ))

        traces.append(go.Scatter3d(
                x=df5[df5["residue_name"] == "HEM"]["x_coordinate"],
                y=df5[df5["residue_name"] == "HEM"]["y_coordinate"],
                z=df5[df5["residue_name"] == "HEM"]["z_coordinate"],
                name="Heme",
                mode="markers",
                marker=dict(
                    color=aa_color_picker("HEM"),
                    line=dict(width=1)
                )
        ))
        return {
            "data": traces,
            "layout": go.Layout(
                height=700,
            )
        }

    ##### ATOMS #####
    elif pdb_value == "v1a3n" and value_type == "atoms":
        traces = []
        for atom in df1[df1["type"] == "ATOM"]["element_symbol"].unique():
            traces.append(go.Scatter3d(
                x=df1[df1["element_symbol"] == atom]["x_coordinate"],
                y=df1[df1["element_symbol"] == atom]["y_coordinate"],
                z=df1[df1["element_symbol"] == atom]["z_coordinate"],
                name=atom,
                mode="markers",
                marker=dict(
                    color=color_picker(atom),
                    line=dict(
                        width=1,
                        color="rgb(0,0,0)",
                    )
                )))
        for atom in df1[df1["type"] == "HETATM"]["element_symbol"].unique():
            traces.append(go.Scatter3d(
                x=df1[(df1["element_symbol"] == atom) & (df1["residue_name"] == "HEM")]["x_coordinate"],
                y=df1[(df1["element_symbol"] == atom) & (df1["residue_name"] == "HEM")]["y_coordinate"],
                z=df1[(df1["element_symbol"] == atom) & (df1["residue_name"] == "HEM")]["z_coordinate"],
                name=atom+" (Heme)",
                mode="markers",
                marker=dict(
                        color=hem_atom_color_picker(atom),
                        line=dict(
                            width=1,
                        )
                    )))
        return {
            "data": traces,
            "layout": go.Layout(
                height=700,
            )
        }
    elif pdb_value == "v2gdm" and value_type == "atoms":
        traces = []
        for atom in df2[df2["type"] == "ATOM"]["element_symbol"].unique():
            traces.append(go.Scatter3d(
                x=df2[df2["element_symbol"] == atom]["x_coordinate"],
                y=df2[df2["element_symbol"] == atom]["y_coordinate"],
                z=df2[df2["element_symbol"] == atom]["z_coordinate"],
                name=atom,
                mode="markers",
                marker=dict(
                    color=color_picker(atom),
                    line=dict(
                        width=1,
                        color="rgb(0,0,0)",
                    )
                )))
        for atom in df2[df2["type"] == "HETATM"]["element_symbol"].unique():
            traces.append(go.Scatter3d(
                x=df2[(df2["element_symbol"] == atom) & (df2["residue_name"] == "HEM")]["x_coordinate"],
                y=df2[(df2["element_symbol"] == atom) & (df2["residue_name"] == "HEM")]["y_coordinate"],
                z=df2[(df2["element_symbol"] == atom) & (df2["residue_name"] == "HEM")]["z_coordinate"],
                name=atom+" (Heme)",
                mode="markers",
                marker=dict(
                        color=hem_atom_color_picker(atom),
                        line=dict(
                            width=1,
                        )
                    )))
        return {
            "data": traces,
            "layout": go.Layout(
                height=700,
            )
        }
    elif pdb_value == "v2qsp" and value_type == "atoms":
        traces = []
        for atom in df3[df3["type"] == "ATOM"]["element_symbol"].unique():
            traces.append(go.Scatter3d(
                x=df3[df3["element_symbol"] == atom]["x_coordinate"],
                y=df3[df3["element_symbol"] == atom]["y_coordinate"],
                z=df3[df3["element_symbol"] == atom]["z_coordinate"],
                name=atom,
                mode="markers",
                marker=dict(
                    color=color_picker(atom),
                    line=dict(
                        width=1,
                        color="rgb(0,0,0)",
                    )
                )))
        for atom in df3[df3["type"] == "HETATM"]["element_symbol"].unique():
            traces.append(go.Scatter3d(
                x=df3[(df3["element_symbol"] == atom) & (df3["residue_name"] == "HEM")]["x_coordinate"],
                y=df3[(df3["element_symbol"] == atom) & (df3["residue_name"] == "HEM")]["y_coordinate"],
                z=df3[(df3["element_symbol"] == atom) & (df3["residue_name"] == "HEM")]["z_coordinate"],
                name=atom+" (Heme)",
                mode="markers",
                marker=dict(
                        color=hem_atom_color_picker(atom),
                        line=dict(
                            width=1,
                        )
                    )))
        return {
            "data": traces,
            "layout": go.Layout(
                height=700,
            )
        }
    elif pdb_value == "v3hrw" and value_type == "atoms":
        traces = []
        for atom in df4[df4["type"] == "ATOM"]["element_symbol"].unique():
            traces.append(go.Scatter3d(
                x=df4[df4["element_symbol"] == atom]["x_coordinate"],
                y=df4[df4["element_symbol"] == atom]["y_coordinate"],
                z=df4[df4["element_symbol"] == atom]["z_coordinate"],
                name=atom,
                mode="markers",
                marker=dict(
                    color=color_picker(atom),
                    line=dict(
                        width=1,
                        color="rgb(0,0,0)",
                    )
                )))
        for atom in df4[df4["type"] == "HETATM"]["element_symbol"].unique():
            traces.append(go.Scatter3d(
                x=df4[(df4["element_symbol"] == atom) & (df4["residue_name"] == "HEM")]["x_coordinate"],
                y=df4[(df4["element_symbol"] == atom) & (df4["residue_name"] == "HEM")]["y_coordinate"],
                z=df4[(df4["element_symbol"] == atom) & (df4["residue_name"] == "HEM")]["z_coordinate"],
                name=atom+" (Heme)",
                mode="markers",
                marker=dict(
                        color=hem_atom_color_picker(atom),
                        line=dict(
                            width=1,
                        )
                    )))
        return {
            "data": traces,
            "layout": go.Layout(
                height=700,
            )
        }
    elif pdb_value == "v4odc" and value_type == "atoms":
        traces = []
        for atom in df5[df5["type"] == "ATOM"]["element_symbol"].unique():
            traces.append(go.Scatter3d(
                x=df5[df5["element_symbol"] == atom]["x_coordinate"],
                y=df5[df5["element_symbol"] == atom]["y_coordinate"],
                z=df5[df5["element_symbol"] == atom]["z_coordinate"],
                name=atom,
                mode="markers",
                marker=dict(
                    color=color_picker(atom),
                    line=dict(
                        width=1,
                        color="rgb(0,0,0)",
                    )
                )))
        for atom in df5[df5["type"] == "HETATM"]["element_symbol"].unique():
            traces.append(go.Scatter3d(
                x=df5[(df5["element_symbol"] == atom) & (df5["residue_name"] == "HEM")]["x_coordinate"],
                y=df5[(df5["element_symbol"] == atom) & (df5["residue_name"] == "HEM")]["y_coordinate"],
                z=df5[(df5["element_symbol"] == atom) & (df5["residue_name"] == "HEM")]["z_coordinate"],
                name=atom+" (Heme)",
                mode="markers",
                marker=dict(
                        color=hem_atom_color_picker(atom),
                        line=dict(
                            width=1,
                        )
                    )))
        return {
            "data": traces,
            "layout": go.Layout(
                height=700,
            )
        }

    ##### SUBUNITS #####
    elif pdb_value == "v1a3n" and value_type == "subunits":
        traces = []
        for subunit in df1[df1["type"] == "ATOM"]["chain_identifier"].unique():
            traces.append(go.Scatter3d(
                x=df1[df1["chain_identifier"] == subunit]["x_coordinate"],
                y=df1[df1["chain_identifier"] == subunit]["y_coordinate"],
                z=df1[df1["chain_identifier"] == subunit]["z_coordinate"],
                name=subunit,
                mode="markers",
                marker=dict(
                    line=dict(width=1)
                )
            ))

        traces.append(go.Scatter3d(
                x=df1[df1["residue_name"] == "HEM"]["x_coordinate"],
                y=df1[df1["residue_name"] == "HEM"]["y_coordinate"],
                z=df1[df1["residue_name"] == "HEM"]["z_coordinate"],
                name="Heme",
                mode="markers",
                marker=dict(
                    color=aa_color_picker("HEM"),
                    line=dict(width=1)
                )
        ))
        return {
            "data": traces,
            "layout": go.Layout(
                height=700,
            )
        }
    elif pdb_value == "v2gdm" and value_type == "subunits":
        traces = []
        for subunit in df2[df2["type"] == "ATOM"]["chain_identifier"].unique():
            traces.append(go.Scatter3d(
                x=df2[df2["chain_identifier"] == subunit]["x_coordinate"],
                y=df2[df2["chain_identifier"] == subunit]["y_coordinate"],
                z=df2[df2["chain_identifier"] == subunit]["z_coordinate"],
                name=subunit,
                mode="markers",
                marker=dict(
                    line=dict(width=1)
                )
            ))

        traces.append(go.Scatter3d(
                x=df2[df2["residue_name"] == "HEM"]["x_coordinate"],
                y=df2[df2["residue_name"] == "HEM"]["y_coordinate"],
                z=df2[df2["residue_name"] == "HEM"]["z_coordinate"],
                name="Heme",
                mode="markers",
                marker=dict(
                    color=aa_color_picker("HEM"),
                    line=dict(width=1)
                )
        ))
        return {
            "data": traces,
            "layout": go.Layout(
                height=700,
            )
        }
    elif pdb_value == "v2qsp" and value_type == "subunits":
        traces = []
        for subunit in df3[df3["type"] == "ATOM"]["chain_identifier"].unique():
            traces.append(go.Scatter3d(
                x=df3[df3["chain_identifier"] == subunit]["x_coordinate"],
                y=df3[df3["chain_identifier"] == subunit]["y_coordinate"],
                z=df3[df3["chain_identifier"] == subunit]["z_coordinate"],
                name=subunit,
                mode="markers",
                marker=dict(
                    line=dict(width=1)
                )
            ))

        traces.append(go.Scatter3d(
                x=df3[df3["residue_name"] == "HEM"]["x_coordinate"],
                y=df3[df3["residue_name"] == "HEM"]["y_coordinate"],
                z=df3[df3["residue_name"] == "HEM"]["z_coordinate"],
                name="Heme",
                mode="markers",
                marker=dict(
                    color=aa_color_picker("HEM"),
                    line=dict(width=1)
                )
        ))
        return {
            "data": traces,
            "layout": go.Layout(
                height=700,
            )
        }
    elif pdb_value == "v3hrw" and value_type == "subunits":
        traces = []
        for subunit in df4[df4["type"] == "ATOM"]["chain_identifier"].unique():
            traces.append(go.Scatter3d(
                x=df4[df4["chain_identifier"] == subunit]["x_coordinate"],
                y=df4[df4["chain_identifier"] == subunit]["y_coordinate"],
                z=df4[df4["chain_identifier"] == subunit]["z_coordinate"],
                name=subunit,
                mode="markers",
                marker=dict(
                    line=dict(width=1)
                )
            ))

        traces.append(go.Scatter3d(
                x=df4[df4["residue_name"] == "HEM"]["x_coordinate"],
                y=df4[df4["residue_name"] == "HEM"]["y_coordinate"],
                z=df4[df4["residue_name"] == "HEM"]["z_coordinate"],
                name="Heme",
                mode="markers",
                marker=dict(
                    color=aa_color_picker("HEM"),
                    line=dict(width=1)
                )
        ))
        return {
            "data": traces,
            "layout": go.Layout(
                height=700,
            )
        }
    elif pdb_value == "v4odc" and value_type == "subunits":
        traces = []
        for subunit in df5[df5["type"] == "ATOM"]["chain_identifier"].unique():
            traces.append(go.Scatter3d(
                x=df5[df5["chain_identifier"] == subunit]["x_coordinate"],
                y=df5[df5["chain_identifier"] == subunit]["y_coordinate"],
                z=df5[df5["chain_identifier"] == subunit]["z_coordinate"],
                name=subunit,
                mode="markers",
                marker=dict(
                    line=dict(width=1)
                )
            ))

        traces.append(go.Scatter3d(
                x=df5[df5["residue_name"] == "HEM"]["x_coordinate"],
                y=df5[df5["residue_name"] == "HEM"]["y_coordinate"],
                z=df5[df5["residue_name"] == "HEM"]["z_coordinate"],
                name="Heme",
                mode="markers",
                marker=dict(
                    color=aa_color_picker("HEM"),
                    line=dict(width=1)
                )
        ))
        return {
            "data": traces,
            "layout": go.Layout(
                height=700,
            )
        }

@app.callback(
    Output("first-species-selected", "children"),
    [Input("species-name", "value")]
)
def set_chosen_name(selected_value):
    for key, value in option_species.items():
        if value == selected_value:
            return dcc.Markdown("You've selected **\"{}\"**".format(key))


@app.callback(
    Output("composition-graph", "figure"),
    [Input("composition-radio", "value"),
     Input("species-name", "value"),
     Input("second-species", "value")]
)
def update_comp_graph(comp_entry, drop_entry, radio_entry):
    name1 = [name for name, value in option_species.items() if value == drop_entry]
    name2 = [name for name, value in option_species.items() if value == radio_entry]
    if comp_entry == "count":
        trace1 = choose_count(drop_entry)
        trace2 = choose_count(radio_entry)
        fig = tools.make_subplots(
            rows=2,
            cols=1,
            shared_xaxes=True,
            subplot_titles=(name1[0], name2[0])
        )
        for trace in trace1:
            fig.append_trace(trace, 1, 1)
        for trace in trace2:
            fig.append_trace(trace, 2, 1)

        fig["layout"].update(
            title="Amino acids count comparison between <br>{} and {}".format(name1[0], name2[0])
        )
        fig["layout"]["yaxis1"].update(title="number of AA")
        fig["layout"]["yaxis2"].update(title="number of AA")
        return fig

    elif comp_entry == "percentage":
        trace1 = choose_percentage(drop_entry)
        trace2 = choose_percentage(radio_entry)
        fig = tools.make_subplots(
            rows=2,
            cols=1,
            shared_xaxes=True,
            subplot_titles=(name1[0], name2[0])
        )
        for trace in trace1:
            fig.append_trace(trace, 1, 1)
        for trace in trace2:
            fig.append_trace(trace, 2, 1)

        fig["layout"].update(
            title="Amino acids percentage comparison between <br>{} and {}".format(name1[0], name2[0])
        )
        fig["layout"]["yaxis1"].update(title="%")
        fig["layout"]["yaxis2"].update(title="%")
        return fig
    elif comp_entry == "weight":
        trace1 = choose_weight(drop_entry)
        trace2 = choose_weight(radio_entry)
        fig = tools.make_subplots(
            rows=2,
            cols=1,
            shared_xaxes=True,
            subplot_titles=(name1[0], name2[0])
        )
        for trace in trace1:
            fig.append_trace(trace, 1, 1)
        for trace in trace2:
            fig.append_trace(trace, 2, 1)

        fig["layout"].update(
            title="Chain weight comparison between <br>{} and {}".format(name1[0], name2[0])
        )
        fig["layout"]["yaxis1"].update(title="g/mol")
        fig["layout"]["yaxis2"].update(title="g/mol")
        return fig


if __name__ == '__main__':
    app.run_server(debug=True)