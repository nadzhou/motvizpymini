import base64
import datetime
import io

import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import dash_table
from stats import seq_extract
from stats import Analysis

import pandas as pd
import numpy as np
import time

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
server = app.server
app.layout = html.Div(
    children=[
        html.Meta(name='viewport', content='width=device-width, initial-scale=1.0'),
        
        html.Div([
                html.Img(
                    src="https://image.freepik.com/free-vector/abstract-dna-symbol_23-2147504932.jpg",
                    className='two columns',
                    style={
                        'height': '`7%',
                        'width': '7%',
                        'float': 'right',
                        'position': 'relative',
                        'margin-top': 0,
                    },
                ),
            html.H1('MotvizPy Mini - Conservation Score Calculator', 
                    style={'textAlign' : 'center', 
                           'font-family' : 'Courier New, monospace'
                           }, 
                    className="inline columns"),
            
            html.H4("This tool uses Shannon Entropy[1] per vertical position to calculate conservation", 
                    style={'textAlign' : 'center', 
                           'font-family' : 'Courier New, monospace'},
                    className="inline columns"
                    ),

            
            html.Div([
                dcc.Upload(
                    id='upload-data',
                    children=html.Div([
                        'Drag and Drop or ',
                        html.Button('Select Alignment File(S)', 
                                    style={
                                        'font-family' : 'Courier New, monospace',
                                        'color' : 'blue'
                                           })
                    ]),style={
                        'width': '100%',
                        'height': '60px',
                        'lineHeight': '60px',
                        'borderWidth': '1px',
                        'borderStyle': 'dashed',
                        'borderRadius': '5px',
                        'textAlign': 'center',
                        'padding': '10px', 
                        'font-family' : 'Courier New, monospace'
                    },
                    # Allow multiple files to be uploaded
                    multiple=True, 
                    loading_state = dict(
                        is_loading=True
                    )
                ),
                html.Div(id='output-data-upload'),
                
            ])
        ]),
    #dcc.Loaddcc.Loading(id="loading-1", children=[html.Div(id="loading-output-1")], type="circle"),
        html.Div(children=[
            html.Footer('1. Batista et al, "An entropy-based approach for the identification of phylogenetically informative genomic regions of Papillomavirus"', 
                    style={'textAlign' : 'left', 
                           'font-family' : 'Courier New, monospace', 
                           'fontSize' : 13}), 
            dcc.Link('Link to paper', href='https://www.sciencedirect.com/science/article/pii/S1567134811003236',
                     style={'textAlign' : 'left', 
                           'font-family' : 'Courier New, monospace', 
                           'fontSize' : 13})
        ])
    ]
) 

def check_file(filename): 
    print(filename)
    supported_files = ["aln", "sth", "phylip", "clustal"]
    
    for ext in supported_files: 
        if ext in filename: 
            return True, ext
    return False, 0

def parse_contents(contents, filename, date):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    
    checked_files, ext = check_file(filename)
    
    try: 
        if checked_files: 
            # Assume that the user uploaded a CSV file
            seq = seq_extract(io.StringIO(decoded.decode('utf-8')), ext)
            seq = [[x for x in y] for y in seq]
            c = Analysis(seq, "1xef")
            c_ent = c.conservation_score(c.seq2np())
            norm_data = c.normalize_data(c_ent)
            norm_data_len = np.arange(len(norm_data))
            
        return html.Div([
            html.Div(children=[
                html.H1(children='Conservation score per nucleotide/amino acid position', 
                        style={
                            'textAlign': 'center',
                            'font-family' : 'Courier New, monospace'
                            }
                        ),

                
                html.Div(children='Interactive representation of conservation score.', 
                    style={
                            'textAlign': 'center',
                            'font-family' : 'Courier New, monospace'
                            }     
                    ),

                dcc.Graph(
                    id='cons-graph',
                    animate=True, 
                    figure={
                        'data': [
                            {'x': norm_data_len, 'y': norm_data, 'type': 'line', 'name': 'SF'},
                        ],
                        'layout': dict(
                            xaxis={'title': 'Sequence Position', 
    },
                            yaxis={'title': 'Conservation score'},
                            hovermode='Closest'
                        )
                    },
                )
            ], className="six columns")
        ],className='row'
    )
                    
    except Exception as e:
        print(e)
        return html.Div(children=
        [
            html.H2(children=
            'ERROR: File upload failed. \
            Please upload file with ".aln" (clustalx files), \
            ".sth" (stockholm) or ".philip" extension',
            style={
                    'textAlign': 'center',
                    'font-family' : 'Courier New, monospace', 
                    'color' : 'red'
                }
            )
        ])



@app.callback(Output('output-data-upload', 'children'),
              [Input('upload-data', 'contents')],
              [State('upload-data', 'filename'),
               State('upload-data', 'last_modified')])


def update_output(list_of_contents, list_of_names, list_of_dates):
    if list_of_contents is not None:
        children = [
            parse_contents(c, n, d) for c, n, d in
            zip(list_of_contents, list_of_names, list_of_dates)]
        return children



    

if __name__ == '__main__':
    app.run_server(debug=True)
