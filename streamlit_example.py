import plotly.express as px
import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots

st.title('Drug Synergy')

df = pd.read_csv('data/average_dose.csv')
compiled = pd.read_csv('data/compiled_plates.csv')
name1 = 'Tazemetostat_1'
name2 = 'JQ1_2'
experiment = 'EOL1_ASSAY_ID_10946'
name1 = st.selectbox('Choose Drug 1', set(df['Name1']))

data1 = df[df['Name1']== name1]
name2 = st.selectbox('Choose Drug 2', set(data1['Name2']))
st.write('Drug1 you selected is:', name1)
st.write('Drug2 you selected is:', name2)
data2 = data1[data1['Name2'] == name2]
experiment = st.selectbox('Choose Experiment', set(data2['Experiment']))
st.write('Experiment you selected is:', experiment)
df1 =data2[data2['Experiment'] == experiment]
df1 = df1[['Name', 'Name1', 'Name2', 'Drug', 'Average']]

#dose response plot
#prediction
drug1_min = df1[df1['Name'] == name1].Drug.min()
drug1_max = df1[df1['Name'] == name1].Drug.max()
glm_model1 = smf.glm('Average ~ Drug', df1[df1['Name'] == name1], family=sm.families.Binomial()).fit()

X1 = np.linspace(drug1_min,drug1_max,100)
Y1 = glm_model1.predict(exog=dict(Drug=X1))

drug2_min = df1[df1['Name'] == name2].Drug.min()
drug2_max = df1[df1['Name'] == name2].Drug.max()
glm_model2 = smf.glm('Average ~ Drug', df1[df1['Name'] == name2], family=sm.families.Binomial()).fit()
X2 = np.linspace(drug2_min,drug2_max,100)
Y2 = glm_model2.predict(exog=dict(Drug=X2))


dose_response = make_subplots(rows=1, cols = 2,
                    x_title='Concentration (uM)',
                    y_title='Percent Survived',
                    vertical_spacing=0.1,
                    horizontal_spacing=0.085,
                    subplot_titles = [name2, name1]
                    )

dose_response.add_trace(go.Scatter(x=df1[df1['Name'] == name1]['Drug'],
                         y = df1[df1['Name'] == name1]['Average'],
                         mode = 'markers',
                         marker_color='red',
                         showlegend = False),
              row=1, col=1)
dose_response.add_trace(go.Scatter(x=X1,
                         y=Y1,
                         marker_color='red',
                         name=name2),
              row=1, col=1)

dose_response.add_trace(go.Scatter(x=df1[df1['Name'] == name2]['Drug'],
                         y = df1[df1['Name'] == name2]['Average'],
                         mode = 'markers',
                         marker_color='blue',
                         showlegend = False),
              row=1, col=2)
dose_response.add_trace(go.Scatter(x=X2,
                         y=Y2,
                         marker_color='blue',
                         name=name1),
              row=1, col=2)

dose_response.update_traces(marker=dict(size=8,
                              line=dict(width=2,
                                        color='DarkSlateGrey')),
                  selector=dict(mode='markers'))
#change axes ranges to be consistent
xmin = min([drug1_min,drug2_min])
xmax = max([drug1_max,drug2_max])
dose_response.update_xaxes(range=[xmin-(xmax-xmin)/20,xmax+(xmax-xmin)/20],tick0 = xmin, tickfont_size = 16)
dose_response.update_yaxes(range=[-0.05,1.05], tick0=0, tickformat = ',.0%', tickfont_size = 16)
dose_response.update_layout(legend=dict(font_size = 20, title_text = 'Treatment'),
                title = dict(text='Individual Dose Response Curves', x=0.5, font_size=30),
                legend_tracegroupgap= 30,
                showlegend=False)
#change axis title size
dose_response.layout.annotations[0]["font"] = {'size': 20}
dose_response.layout.annotations[1]["font"] = {'size': 20}
#dose_response.show()
st.plotly_chart(dose_response)

# ##############fig 2################
# #IC50 plots
compiled_subset = compiled[(compiled['Name1'] == name1) & (compiled['Name2'] == name2) & (compiled['Experiment'] == experiment)]
compiled_avg = compiled_subset.groupby(['Name1','Drug1','Group1','Name2','Drug2','Group2'], as_index=False).agg({'Value': ['mean', 'std']})
compiled_avg.columns = ['Name1','Drug1','Group1','Name2','Drug2','Group2', 'Mean', 'sd']

#construct color scheme
color_palette=['rgb(227,74,51)']
blues =  compiled_avg[['Group1', 'Group2']].max().max()-1
start_blue = [int(x) for x in px.colors.sequential.Blues[2].replace('rgb(','').replace(")","").split(",")]
end_blue = [int(x) for x in px.colors.sequential.Blues[-1].replace('rgb(','').replace(")","").split(",")]
converted_color = np.vstack([ np.linspace(x[0], x[1], blues) for x in zip(start_blue, end_blue)]).astype(int)
for i in range(0, converted_color.shape[1]):
    color_palette.append("rgb(" + ",".join([str(x) for x in converted_color[:,i]]) + ")")

#Drug1:
drug1_min = compiled_avg.Drug1.min()
drug1_max = compiled_avg.Drug1.max()
X1 = np.linspace(drug1_min,drug1_max,100)
Y1 = np.zeros(shape=(len(X1), len(set(compiled_avg.Group2))))

for i in set(compiled_avg.Group2):
    temp = compiled_avg[compiled_avg.Group2 == i]
    glm_model = smf.glm('Mean ~ Drug1', temp, family=sm.families.Binomial()).fit()
    y = glm_model.predict(exog=dict(Drug1=X1))
    y[(X1 < temp.Drug1.min()) & (X1 > temp.Drug1.max())] = np.nan
    Y1[:,i-1] = np.array(y)

#Drug2:
drug2_min = compiled_avg.Drug2.min()
drug2_max = compiled_avg.Drug2.max()
X2 = np.linspace(drug2_min,drug2_max,100)
Y2 = np.zeros(shape=(len(X2), len(set(compiled_avg.Group1))))

for i in set(compiled_avg.Group1):
    temp = compiled_avg[compiled_avg.Group1 == i]
    glm_model = smf.glm('Mean ~ Drug2', temp, family=sm.families.Binomial()).fit()
    y = glm_model.predict(exog=dict(Drug2=X2))
    y[(X2 < temp.Drug2.min()) & (X2 > temp.Drug2.max())] = np.nan
    Y2[:,i-1] = np.array(y)
#xaxis range
xmin = min([drug1_min,drug2_min])
xmax = max([drug1_max,drug2_max])
#figure
IC50 = make_subplots(rows=1, cols = 2,
                    x_title='Concentration (uM)',
                    y_title='Percent Survived',
                    vertical_spacing=0.1,
                    horizontal_spacing=0.085,
                    subplot_titles = [name1, name2]
                    )
#scatter plots
IC50.add_trace(go.Scatter(x=compiled_avg['Drug1'],
                        y = compiled_avg['Mean'],
                        mode = 'markers',
                        marker_color=compiled_avg['Group2'],
                        marker_colorscale=color_palette,
                        showlegend = False),
              row=1, col=1)

IC50.add_trace(go.Scatter(x=compiled_avg['Drug2'],
                        y = compiled_avg['Mean'],
                        mode = 'markers',
                        marker_color=compiled_avg['Group1'],
                        marker_colorscale=color_palette,
                        showlegend = False),
              row=1, col=2)

#logistic regression line
for j in range(0,Y1.shape[1]):
    IC50.add_trace(go.Scatter(x=X1,
                            y= Y1[:,j],
                            mode = "lines",
                            line_color = color_palette[j],
                            showlegend = False),
                  row=1,col=1)

for j in range(0,Y2.shape[1]):
    IC50.add_trace(go.Scatter(x=X2,
                            y= Y2[:,j],
                            mode = "lines",
                            line_color = color_palette[j],
                            name = j+1),
                  row=1,col=2)

#adjustments
#IC50.update_xaxes(range=[xmin-(xmax-xmin)/20,xmax+(xmax-xmin)/20],tick0 = xmin, tickfont_size = 16)
IC50.update_yaxes(range=[-0.05,1.05], tick0=0, tickformat = ',.0%', tickfont_size = 16)
IC50.update_layout(legend=dict(font_size = 20, title_text='Dose Combination'),
                title = dict(text='Shifting IC50 Curves', x=0.5, font_size=30),
                legend_tracegroupgap= 30)
#change axis title size
IC50.layout.annotations[0]["font"] = {'size': 20}
IC50.layout.annotations[1]["font"] = {'size': 20}
#IC50.show()
st.plotly_chart(IC50)