# import dash
# import dash_core_components as dcc
# import dash_html_components as html
# import plotly
# from plotly.offline import iplot, init_notebook_mode
# init_notebook_mode(connected = True)

# external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

# app = dash.Dash(__name__, external_stylesheets=external_stylesheets)


# gene='CTCF'
# chr='chr1'
# cell='A-549'



# app.layout = html.Div(children=[
#     html.H1(children='Hello Dash'),

#     html.Div([
#     html.Label('Cell Line'),
#     dcc.Input(value='A-549', type='text'),
#     html.Label('Chromosome'),
#     dcc.Input(value='chr1', type='text'),
#     html.Label('Transcription Factor'),
#     dcc.Input(value='CTCF', type='text'),
#     ],style={'width': '49%', 'display': 'inline-block'}),

#     indices = [i for i, s in enumerate(traces) if cell+'_'+gene in s]
#     data=pd.read_csv(traces[indices[0]],sep='\t',usecols=[0,1,2,3,4,5,6,7,8,12,13],names=["chr", "start", "end",'weight',"wgbs",'gene',"ChIPTF",'shW','shWG','array','shA']) 
#     data=data[data['weight']!=data['wgbs']] ##subset of motif or entire motif flag
#     data=data[data['weight']!=data['array']]
#     data['cell']=os.path.basename(traces[indices[0]]).split('_')[0]
#     data=data[data['chr']==chr]

#     dcc.Graph(
#         id='example-graph',
#         figure={
#             'data': [
#                 dict(
#                     x= np.log10(data.start),#df[df['continent'] == i]['gdp per capita'],
#                     y= data.weight,#df[df['continent'] == i]['life expectancy'],
#                     text=data.cell+'PWM='+data.weight, #df[df['continent'] == i]['country'],
#                     mode='markers',
#                     opacity=0.7,
#                     marker={
#                         'size': 15,
#                         'line': {'width': 0.5, 'color': 'white'}
#                     },
#                     name=i
#                 ) for i in df.continent.unique()
#             ],
#             'layout': dict(
#                 xaxis={'type': 'log', 'title': 'Location on Chromosome'},
#                 yaxis={'title': 'PWM weight'},
#                 margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
#                 legend={'x': 0, 'y': 1},
#                 hovermode='closest'
#             )
#         }


# # fig = go.Figure()#make_subplots(rows=4, cols=1)
# # x=np.log10(data.start)

# # fig.add_trace(go.Bar(x=np.log10(data.start),y=data.weight))#, row=1, col=1)
# # fig.add_trace(go.Bar(x=np.log10(data.start),y=1-data.wgbs))#, row=1, col=1)
# # fig.add_trace(go.Bar(x=np.log10(data.start),y=1-data.array))#, row=1, col=1)
# # fig.add_trace(go.Bar(x=np.log10(data.start),y=data.ChIPTF))#, row=1, col=1)
# # fig.update_layout(barmode='relative')#,height=600, width=600, title_text="Stacked Subplots")
# # fig.show()
        

#     )
# ])

# if __name__ == '__main__':
#     app.run_server(debug=True)




import dash
import glob
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
traces= glob.glob('data/MotifPipeline/MeArrayIntersectFULLshuf/*')
indices = [i for i, s in enumerate(traces) if cell+'_'+gene in s]
df=pd.read_csv(traces[indices[0]],sep='\t',usecols=[0,1,2,3,4,5,6,7,8,12,13],names=["chr", "start", "end",'weight',"wgbs",'gene',"ChIPTF",'shW','shWG','array','shA']) 
df=df[df['weight']!=df['wgbs']] ##subset of motif or entire motif flag
df=df[df['weight']!=df['array']]
df['cell']=os.path.basename(traces[indices[0]]).split('_')[0]
# df = #pd.read_csv('https://plotly.github.io/datasets/country_indicators.csv')

available_indicators = df['gene'].unique()

app.layout = html.Div([
    html.Div([

        html.Div([
            dcc.Dropdown(
                id='crossfilter-xaxis-column',
                options=[{'label': i, 'value': i} for i in available_indicators],
                value='Fertility rate, total (births per woman)'
            ),
            dcc.RadioItems(
                id='crossfilter-xaxis-type',
                options=[{'label': i, 'value': i} for i in ['Linear', 'Log']],
                value='Linear',
                labelStyle={'display': 'inline-block'}
            )
        ],
        style={'width': '49%', 'display': 'inline-block'}),

        html.Div([
            dcc.Dropdown(
                id='crossfilter-yaxis-column',
                options=[{'label': i, 'value': i} for i in available_indicators],
                value='Life expectancy at birth, total (years)'
            ),
            dcc.RadioItems(
                id='crossfilter-yaxis-type',
                options=[{'label': i, 'value': i} for i in ['Linear', 'Log']],
                value='Linear',
                labelStyle={'display': 'inline-block'}
            )
        ], style={'width': '49%', 'float': 'right', 'display': 'inline-block'})
    ], style={
        'borderBottom': 'thin lightgrey solid',
        'backgroundColor': 'rgb(250, 250, 250)',
        'padding': '10px 5px'
    }),

    html.Div([
        dcc.Graph(
            id='crossfilter-indicator-scatter',
            hoverData={'points': [{'customdata': 'Japan'}]}
        )
    ], style={'width': '49%', 'display': 'inline-block', 'padding': '0 20'}),
    html.Div([
        dcc.Graph(id='x-time-series'),
        dcc.Graph(id='y-time-series'),
    ], style={'display': 'inline-block', 'width': '49%'}),

    html.Div(dcc.Slider(
        id='crossfilter-cell--slider',
        min=df['Year'].min(),
        max=df['Year'].max(),
        value=df['Year'].max(),
        marks={str(cell): str(cell) for cell in df['Cell'].unique()},
        step=None
    ), style={'width': '49%', 'padding': '0px 20px 20px 20px'})
])


@app.callback(
    dash.dependencies.Output('crossfilter-indicator-scatter', 'figure'),
    [dash.dependencies.Input('crossfilter-xaxis-column', 'value'),
     dash.dependencies.Input('crossfilter-yaxis-column', 'value'),
     dash.dependencies.Input('crossfilter-xaxis-type', 'value'),
     dash.dependencies.Input('crossfilter-yaxis-type', 'value'),
     dash.dependencies.Input('crossfilter-year--slider', 'value')])
def update_graph(xaxis_column_name, yaxis_column_name,
                 xaxis_type, yaxis_type,
                 year_value):
    dff = df[df['Year'] == year_value]

    return {
        'data': [dict(
            x=dff[dff['Indicator Name'] == xaxis_column_name]['Value'],
            y=dff[dff['Indicator Name'] == yaxis_column_name]['Value'],
            text=dff[dff['Indicator Name'] == yaxis_column_name]['Country Name'],
            customdata=dff[dff['Indicator Name'] == yaxis_column_name]['Country Name'],
            mode='markers',
            marker={
                'size': 15,
                'opacity': 0.5,
                'line': {'width': 0.5, 'color': 'white'}
            }
        )],
        'layout': dict(
            xaxis={
                'title': xaxis_column_name,
                'type': 'linear' if xaxis_type == 'Linear' else 'log'
            },
            yaxis={
                'title': yaxis_column_name,
                'type': 'linear' if yaxis_type == 'Linear' else 'log'
            },
            margin={'l': 40, 'b': 30, 't': 10, 'r': 0},
            height=450,
            hovermode='closest'
        )
    }


def create_time_series(dff, axis_type, title):
    return {
        'data': [dict(
            x=dff['Year'],
            y=dff['Value'],
            mode='lines+markers'
        )],
        'layout': {
            'height': 225,
            'margin': {'l': 20, 'b': 30, 'r': 10, 't': 10},
            'annotations': [{
                'x': 0, 'y': 0.85, 'xanchor': 'left', 'yanchor': 'bottom',
                'xref': 'paper', 'yref': 'paper', 'showarrow': False,
                'align': 'left', 'bgcolor': 'rgba(255, 255, 255, 0.5)',
                'text': title
            }],
            'yaxis': {'type': 'linear' if axis_type == 'Linear' else 'log'},
            'xaxis': {'showgrid': False}
        }
    }


@app.callback(
    dash.dependencies.Output('x-time-series', 'figure'),
    [dash.dependencies.Input('crossfilter-indicator-scatter', 'hoverData'),
     dash.dependencies.Input('crossfilter-xaxis-column', 'value'),
     dash.dependencies.Input('crossfilter-xaxis-type', 'value')])
def update_y_timeseries(hoverData, xaxis_column_name, axis_type):
    country_name = hoverData['points'][0]['customdata']
    dff = df[df['Country Name'] == country_name]
    dff = dff[dff['Indicator Name'] == xaxis_column_name]
    title = '<b>{}</b><br>{}'.format(country_name, xaxis_column_name)
    return create_time_series(dff, axis_type, title)


@app.callback(
    dash.dependencies.Output('y-time-series', 'figure'),
    [dash.dependencies.Input('crossfilter-indicator-scatter', 'hoverData'),
     dash.dependencies.Input('crossfilter-yaxis-column', 'value'),
     dash.dependencies.Input('crossfilter-yaxis-type', 'value')])
def update_x_timeseries(hoverData, yaxis_column_name, axis_type):
    dff = df[df['Country Name'] == hoverData['points'][0]['customdata']]
    dff = dff[dff['Indicator Name'] == yaxis_column_name]
    return create_time_series(dff, axis_type, yaxis_column_name)


if __name__ == '__main__':
    app.run_server(debug=True)