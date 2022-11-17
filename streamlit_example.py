from curses.ascii import DC1
import plotly.express as px
import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots

st.set_page_config(layout="wide")
st.title('Drug Synergy')

df = pd.read_csv('data/average_dose.csv')
df['Name1'] = [x.replace("_1","") for x in df.Name1]
df['Name2'] = [x.replace("_2","") for x in df.Name2]
compiled = pd.read_csv('data/compiled_plates.csv')
compiled['Name1'] = [x.replace("_1","") for x in compiled.Name1]
compiled['Name2'] = [x.replace("_2","") for x in compiled.Name2]

individual_doses = pd.read_csv('data/individual_doses.csv')
individual_doses['Name1'] = [x.replace("_1","") for x in individual_doses.Name1]
individual_doses['Name2'] = [x.replace("_2","") for x in individual_doses.Name2]


name1 = "ORY-1001"
name2 = 'Methotrexate'

experiment = 'MOLM14_ASSAY_ID_8872'
data1 = df[df['Name1']== name1]
data2 = data1[data1['Name2'] == name2]
df1 = data2[data2['Experiment'] == experiment]
df1 = df1[['Name', 'Name1', 'Name2', 'Drug', 'Average']]


name1 = st.selectbox('Choose Drug 1', set(df['Name1']))
data1 = df[df['Name1']== name1]
name2 = st.selectbox('Choose Drug 2', set(data1['Name2']))
st.write('Drug1 you selected is:', name1)
st.write('Drug2 you selected is:', name2)
data2 = data1[data1['Name2'] == name2]
experiment = st.selectbox('Choose Experiment', set(data2['Experiment']))
st.write('Experiment you selected is:', experiment)
df1 = data2[data2['Experiment'] == experiment]
df1 = df1[['Name', 'Name1', 'Name2', 'Drug', 'Average']]

#dose response plot
#prediction
#Name1 = name1
#Name2 = name2
Name1 = name1 + '_1'
Name2 = name2 + '_2'
drug1_min = df1[df1['Name'] == Name1].Drug.min()
drug1_max = df1[df1['Name'] == Name1].Drug.max()
glm_model1 = smf.glm('Average ~ Drug', df1[df1['Name'] == Name1], family=sm.families.Binomial()).fit()

X1 = np.linspace(drug1_min,drug1_max,100)
Y1 = glm_model1.predict(exog=dict(Drug=X1))

drug2_min = df1[df1['Name'] == Name2].Drug.min()
drug2_max = df1[df1['Name'] == Name2].Drug.max()
glm_model2 = smf.glm('Average ~ Drug', df1[df1['Name'] == Name2], family=sm.families.Binomial()).fit()
X2 = np.linspace(drug2_min,drug2_max,100)
Y2 = glm_model2.predict(exog=dict(Drug=X2))

#figure 1, dose_response
dose_response = make_subplots(rows=1, cols = 2,
                    x_title='Concentration (uM)',
                    y_title='Percent Survived',
                    vertical_spacing=0.1,
                    horizontal_spacing=0.085,
                    subplot_titles = [name2, name1]
                    )
#scatter
dose_response.add_trace(go.Scatter(x = df1[df1['Name'] == Name1]['Drug'],
                         y = df1[df1['Name'] == Name1]['Average'],
                         mode = 'markers',
                         marker_color='red',
                         showlegend = False,
                         hovertemplate = '%{y: ,.2%}<extra></extra>',
                         ),
              row=1, col=1)
#regerssion line
dose_response.add_trace(go.Scatter(x=X1,
                         y=Y1,
                         marker_color='red',
                         name=name2,
                        hoverinfo='skip'),
              row=1, col=1)

#plot2
#scatter
dose_response.add_trace(go.Scatter(x=df1[df1['Name'] == Name2]['Drug'],
                         y = df1[df1['Name'] == Name2]['Average'],
                         mode = 'markers',
                         marker_color='blue',
                         showlegend = False,
                         hovertemplate = '%{y: ,.2%}<extra></extra>',
                         ),
              row=1, col=2)
#regression line
dose_response.add_trace(go.Scatter(x=X2,
                         y=Y2,
                         marker_color='blue',
                         name=name1,
                         hoverinfo='skip'),
              row=1, col=2)

dose_response.update_traces(marker=dict(size=8,
                              line=dict(width=2,
                                        color='DarkSlateGrey')),
                  selector=dict(mode='markers'))

#change axes ranges to be consistent
xmin = min([drug1_min,drug2_min])
xmax = max([drug1_max,drug2_max])
#dose_response.update_xaxes(range=[xmin-(xmax-xmin)/20,xmax+(xmax-xmin)/20],tick0 = xmin, tickfont_size = 16)
dose_response.update_yaxes(range=[-0.05,1.05], tick0=0, tickformat = ',.0%', tickfont_size = 16)
dose_response.update_layout(legend=dict(font_size = 20, title_text = 'Treatment'),
                title = dict(text='Individual Dose Response Curves', x=0.5, font_size=30),
                legend_tracegroupgap= 30,
                showlegend=False)
#change axis title size
dose_response.layout.annotations[0]["font"] = {'size': 20}
dose_response.layout.annotations[1]["font"] = {'size': 20}
#dose_response.show()

st.plotly_chart(dose_response,use_container_width=True)

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

# st.write('color palette is ', color_palette)
# st.write('blues is ', blues)
# st.write('start blue is', start_blue)
# st.write('end blue is ', end_blue)


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
    if len(set(temp['Mean'])) > 1:
        glm_model = smf.glm('Mean ~ Drug1', temp, family=sm.families.Binomial()).fit()
        y = glm_model.predict(exog=dict(Drug1=X1))
    else:
        y = np.ones(len(X1)) * list(set(temp['Mean']))
    y[(X1 < temp.Drug1.min()) & (X1 > temp.Drug1.max())] = np.nan
    Y1[:,i-1] = np.array(y)

#Drug2:
drug2_min = compiled_avg.Drug2.min()
drug2_max = compiled_avg.Drug2.max()
X2 = np.linspace(drug2_min,drug2_max,100)
Y2 = np.zeros(shape=(len(X2), len(set(compiled_avg.Group1))))

for i in set(compiled_avg.Group1):
    temp = compiled_avg[compiled_avg.Group1 == i]
    if len(set(temp['Mean'])) > 1:
        glm_model = smf.glm('Mean ~ Drug2', temp, family=sm.families.Binomial()).fit()
        y = glm_model.predict(exog=dict(Drug2=X2))
    else:
        y = np.ones(len(X2)) * list(set(temp['Mean']))
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
                    subplot_titles = [name2, name1]
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
                            showlegend = False,
                            hoverinfo='skip'),
                  row=1,col=1)

for j in range(0,Y2.shape[1]):
    IC50.add_trace(go.Scatter(x=X2,
                            y= Y2[:,j],
                            mode = "lines",
                            line_color = color_palette[j],
                            name = j+1,
                            hoverinfo='skip'),
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

st.plotly_chart(IC50, use_container_width=True)

#FA plots
FAdata = compiled_subset.groupby(['Name1','Drug1','Name2','Drug2'], as_index=False).agg({'Fa': ['mean']})
FAdata.columns = ['Name1','Drug1','Name2','Drug2', 'Fa']

#replace min with 0
drug1_min = FAdata['Drug1'].min()
drug2_min = FAdata['Drug2'].min()
FAdata['Drug1'].replace({drug1_min:0}, inplace=True)
FAdata['Drug2'].replace({drug2_min:0}, inplace=True)

FAdata_wide = pd.pivot(FAdata, index='Drug2', columns='Drug1', values='Fa')

# Fraction affected plot
FAplot = go.Figure(data=go.Heatmap(
                    z=FAdata_wide,
                    x=FAdata_wide.columns.map(str),
                    y=FAdata_wide.index.map(str),
                    text = FAdata_wide,
                    texttemplate="%{text:.2~%}",
                    textfont={"size":16},
                    colorscale = 'Reds',
                    hoverinfo='skip'))

FAplot.update_layout(title = dict(text='Percentage Affected', x=0.5, font_size=30),
                   xaxis_title=name1 + ' Concentration (uM)',
                   yaxis_title=name2 + ' Concentration (uM)')

#FAplot.show()

st.plotly_chart(FAplot, use_container_width=True)

#Median Effect plot
# Get the slope and intercept

individual_subset = individual_doses[(individual_doses['Name1'] == name1) & (individual_doses['Name2'] == name2) & (individual_doses['Experiment'] == experiment)]
d1 = compiled_subset[(compiled_subset['Group2'] ==1) & (compiled_subset['Drug1'] != drug1_min)]
d2 = compiled_subset[(compiled_subset['Group1'] ==1) & (compiled_subset['Drug2'] != drug2_min)]

d1_model = smf.ols(formula='np.log10(FaFu) ~ np.log10(Drug1)', data=d1).fit()
d2_model = smf.ols(formula='np.log10(FaFu) ~ np.log10(Drug2)', data=d2).fit()
#coef d1_model.rsquared, d1_model.params
d1_Dm = 10 ** (-d1_model.params[0]/d1_model.params[1])
d2_Dm = 10 ** (-d2_model.params[0]/d2_model.params[1])

#LM fit
idx1 = (individual_subset['Drug'] == 0) & (individual_subset['Name'] == Name1)
individual_subset.loc[idx1, 'Drug'] = individual_subset.loc[idx1,'min1']/2
idx2 = (individual_subset['Drug'] == 0) & (individual_subset['Name'] == Name2)
individual_subset.loc[idx2, 'Drug'] = individual_subset.loc[idx2,'min2']/2
filter = (individual_subset['Drug1']  == individual_subset['min1'].iloc[0]/2) & (individual_subset['Drug2']  == individual_subset['min2'].iloc[0]/2)
individual_subset = individual_subset.loc[~filter,:]

d1_model = smf.ols(formula='np.log10(FaFu) ~ np.log10(Drug1)', data=individual_subset[individual_subset['Name'] == Name1]).fit()
d2_model = smf.ols(formula='np.log10(FaFu) ~ np.log10(Drug2)', data=individual_subset[individual_subset['Name'] == Name2]).fit()

drug1_min = individual_subset[individual_subset['Name'] == Name1]['Drug'].min()
drug1_max = individual_subset[individual_subset['Name'] == Name1]['Drug'].max()
drug2_min = individual_subset[individual_subset['Name'] == Name2]['Drug'].min()
drug2_max = individual_subset[individual_subset['Name'] == Name2]['Drug'].max()
X1 = np.linspace(drug1_min, drug1_max, 100)
X2 = np.linspace(drug2_min, drug2_max,100)
Y1 = d1_model.predict(exog=dict(Drug1=X1))
Y2 = d2_model.predict(exog=dict(Drug2=X2))

#the figure
med_effect = make_subplots(rows=1, cols = 2,
                    x_title='log10[Dose(uM)]',
                    y_title='log10(Fa/Fu)',
                    vertical_spacing=0.1,
                    horizontal_spacing=0.085,
                    subplot_titles = [name1, name2]
                    )
med_effect.add_trace(go.Scatter(x=np.log10(individual_subset[individual_subset['Name'] == Name1]['Drug']),
                         y = np.log10(individual_subset[individual_subset['Name'] == Name1]['FaFu']),
                         mode = 'markers',
                         marker_color='red',
                         showlegend = False),
              row=1, col=1)
med_effect.add_trace(go.Scatter(x=np.log10(individual_subset[individual_subset['Name'] == Name2]['Drug']),
                         y = np.log10(individual_subset[individual_subset['Name'] == Name2]['FaFu']),
                         mode = 'markers',
                         marker_color='blue',
                         showlegend = False),
              row=1, col=2)
med_effect.add_trace(go.Scatter(x=np.log10(X1),
                         y=Y1,
                         marker_color='red',
                         name=name2),
              row=1, col=1)
med_effect.add_trace(go.Scatter(x=np.log10(X2),
                         y=Y2,
                         marker_color='blue',
                         name=name1),
              row=1, col=2)
med_effect.update_traces(marker=dict(size=8,
                              line=dict(width=2,
                                        color='DarkSlateGrey')),
                  selector=dict(mode='markers'))

#change axes ranges to be consistent
ymin = np.log10(min(individual_subset.FaFu)) - 0.1
ymax = np.log10(max(individual_subset.FaFu)) + 0.1
xmin = np.log10(min(drug1_min,drug2_min))

med_effect.add_trace(go.Scatter(x= [xmin+0.2],
                                y= [ymax-0.2],
                                mode='text',
                                text= 'adj. R<sup>2</sup> = '+ str(np.round(d1_model.rsquared_adj,2)),
                                textfont = dict(color='red',
                                                size=20)),
                     row=1,col=1)

med_effect.add_trace(go.Scatter(x= [xmin+0.2],
                                y= [ymax-0.2],
                                mode='text',
                                text= 'adj. R<sup>2</sup> = '+ str(np.round(d2_model.rsquared_adj,2)),
                                textfont = dict(color='blue',
                                                size=20)),
                     row=1,col=2)

med_effect.update_yaxes(range=[ymin,ymax], tick0=ymin, tickfont_size = 16)
med_effect.update_layout(showlegend=False)
#change axis title size
med_effect.layout.annotations[0]["font"] = {'size': 20}
med_effect.layout.annotations[1]["font"] = {'size': 20}
#med_effect.show()
st.plotly_chart(med_effect, use_container_width=True)


#Isobologram
cell_dat = compiled_subset[(compiled_subset.Drug1 != compiled_subset.min1/2) & (compiled_subset.Drug2 != compiled_subset.min2/2)].groupby(['Drug1', 'Drug2', 'Name1', 'Name2'],as_index=False).agg({'Fa': 'mean', 'Fu': 'mean'})
cell_dat.columns = ['Drug1','Drug2', 'Name1', 'Name2', 'Fa', 'Fu']

# Calculate Dx1 and Dx2
# Dx = Dm[fa/fu]^1/m
cell_dat['Dx1'] = d1_Dm*(cell_dat.Fa/cell_dat.Fu)**(1/d1_model.params[1])
cell_dat['Dx2'] = d2_Dm*(cell_dat.Fa/cell_dat.Fu)**(1/d2_model.params[1])
cell_dat['I1'] = cell_dat.Drug1 / cell_dat.Dx1
cell_dat['I2'] = cell_dat.Drug2 / cell_dat.Dx2
cell_dat['CI'] = cell_dat.I1 + cell_dat.I2

isobol_dat = cell_dat.dropna(subset='CI')

# Range of Drug Combination Indexs ----------------------------------------
# From Figure 4 "Theoretical Basis, Experimental Design, and Computerized Simulation of Synergism and Antagonism in Drug Combination Studies"0.
# pick the last syn_val that CI >, and assign description
syn_val = [0.3, 0.7, 1, 1.45, 3]
syn_desc = {"Strong synergism": 0.3,
            "Synergism": 0.7,
            "Additive": 1,
            "Antagonism": 1.45,
            "Strong antagonism": 3}

isobol_dat.loc[isobol_dat.CI>3, 'CI'] = 3

temp = list(syn_desc.values())
desc = []
for CI in isobol_dat.CI:
    if CI < list(syn_desc.values())[0]:
        d = list(syn_desc.keys())[0]
    elif CI > list(syn_desc.values())[-1]:
        d = list(syn_desc.keys())[-1]
    else:
        id = [idx for idx, val in enumerate(temp) if CI <= val][0]-1
    d = list(syn_desc.keys())[id]
    desc.append(d)

isobol_dat['CI_desc'] = desc

# ci_max = isobol_dat.CI.max() if isobol_dat.CI.max() > 2 else 2
# scale_val = [0,1/ci_max, 1]
mean_CI = np.around(isobol_dat.CI.mean(), decimals= 3)
mean_CI_se = np.around(isobol_dat.CI.std()/np.sqrt(isobol_dat.shape[0]), decimals=3)

#make triangle?
triangle_x = [0,0,1,0]
triangle_y = [1,0,0,1]
#triangle = data.frame(x = c(0,1,0), y = c(0,0,1))

#if isobol_dat.shape()[0] >0
isobol_wide = pd.pivot(isobol_dat, index='Drug2', columns='Drug1', values='CI').dropna()

isobologram = make_subplots(rows=1, cols = 2,
                            subplot_titles=['Normalized Isobologram', 'Range of Combination Indices'])

#Isobologram
#triangle
max_scale = isobol_dat[['I1','I2']].max().max()
if max_scale <=1:
    max_scale = 1.2
elif max_scale <=2:
    max_scale=2.2
else:
    max_scale = 3.5

isobologram.add_trace(go.Scatter(x=triangle_x,
                                y=triangle_y,
                                mode='lines',
                                line_color='black',
                                line_width=3,
                                hoverinfo='skip',showlegend=False),
                      1,1)
#points
isobologram.add_trace(go.Scatter(x=isobol_dat.I1,
                                y = isobol_dat.I2,
                                mode = 'markers',
                                marker = dict(color=np.log10(isobol_dat.CI),
                                            size=20,
                                            line=dict(width=2,
                                                        color='black'),
                                            colorscale='RdBu',
                                            cmin=np.log10(0.3),
                                            cmax=np.log10(3),
                                            cmid=0),
                                customdata=isobol_dat[['CI','CI_desc']],
                                hovertemplate='CI = %{customdata[0]:.2} <br>Category = %{customdata[1]} <extra></extra>',
                                showlegend=False),
                      1,1)
#texts
isobologram.add_trace(go.Scatter(x=[2.7],
                                 y=[3],
                                 mode='text',
                                 text=[ str(mean_CI) + "+/-" + str(mean_CI_se)],
                                 textposition='bottom center',
                                 textfont=dict(family='sans serif',
                                               size=20),
                                 hoverinfo='skip',
                                 showlegend=False),
                      1,1)

isobologram.update_xaxes(range=[-0.2,max_scale],tick0 = 0, dtick=1, tickfont_size = 16, title= 'Compound 1',showline = False, tickwidth=3, ticklen=6, row=1, col=1)
isobologram.update_yaxes(range=[-0.2,max_scale],tick0 = 0, dtick=1, tickfont_size = 16, title= 'Compound 2',showline = False, tickwidth=3, ticklen=6, row=1, col=1)

#CI plot
isobologram.add_trace(go.Heatmap(z=np.log10(isobol_wide),
                                x=isobol_wide.columns.map(str),
                                y=isobol_wide.index.map(str),
                                zmax=np.log10(3), zmin=np.log10(0.3), zmid = 0,
                                text = isobol_wide,
                                texttemplate="%{text:.2}",
                                textfont={"size":16},
                                colorscale='RdBu',
                                hoverinfo='none',
                                colorbar=dict(tickmode='array',
                                            tickvals=np.log10(syn_val),
                                            tickcolor='black',
                                            ticktext=list(syn_desc.keys()))),
                      1,2)

isobologram.update_xaxes(title= 'Compound 2 Concentration (uM)', row=1, col=2)
isobologram.update_yaxes(title= 'Compound 1 Concentration (uM)', row=1, col=2)

isobologram.update_layout(title = dict(x=0.5, font_size=30),
                          template='simple_white')

#isobologram.show()
st.plotly_chart(isobologram, use_container_width=True)