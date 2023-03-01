#!/usr/bin/python

import pandas as pd
import matplotlib.pyplot as plt
import sys
import numpy as np
import scipy as sc
import sklearn
from sklearn import decomposition
from sklearn import cluster
import seaborn


# cells excluded because of low read count
cell_sets = {}
with open(sys.argv[3]) as f:
	for line in f:
		x = line.rstrip().split("\t")
		cell_sets[x[0]] = x[1:]

# get parameters
expression_file = sys.argv[1]
cell_group_file = sys.argv[2]

# set colors for different days and treatments
day_color       = pd.Series(['#fdcc8a','#fc8d59','#e34a33','#b30000'], index=["Day_4", "Day_6", "Day_8", "Day_12"])
day_size        = pd.Series([7,9,11,13], index=["Day_4", "Day_6", "Day_8", "Day_12"])
#day_size = day_size*4
#treatment_color = pd.Series(['#e41a1c','#377eb8','#4daf4a']          , index=["sh733", "sh737", "shCtrl"]) # red, blue, green
treatment_color = pd.Series(['#e41a1c','orange','gray']          , index=["sh733", "sh737", "shCtrl"]) # red, blue, green
cell_groups = pd.read_csv(cell_group_file, sep=",", index_col=0)
cell_groups["day_color"] = day_color[cell_groups["Day"]].values
cell_groups["day_size"] = day_size[cell_groups["Day"]].values
cell_groups["treatment_color"] = treatment_color[cell_groups["Treatment"]].values

# read expression
expression_table = pd.read_csv(expression_file, sep=",").transpose()
# remove low read count cells
expression_table = expression_table.drop(cell_sets["low_read_count_cells"])
# log transform
expression_table += 1
expression_table = expression_table.apply(np.log2)
# remove genes with <10 cells expressing them
min_expression = 0.1
min_cells = 10
expressed_genes = (expression_table > min_expression).sum() > min_cells
expression_table = expression_table.loc[ : , expressed_genes]

# create annotation dataframe
# !!! order of rows MUST correspond to rows in the expression table after all manipulations !!!!!!!!
annotations = pd.DataFrame(index=expression_table.index)
annotations["quality_color"] = "blue"
#annotations.loc[low_read_count_cells, "quality_color"]="red"

annotations["day_color"] = cell_groups.loc[annotations.index, "day_color"]
annotations["treatment_color"] = cell_groups.loc[annotations.index, "treatment_color"]
annotations["day_size"] = cell_groups.loc[annotations.index, "day_size"]
#annotations.loc[lhcb, "color"]="orange"
#annotations.loc[shcb, "color"]="pink"

def annotate_df(row):
	dist = (tsne_transformed_expression_2d - row).abs().sum(axis=1).sort_values()[1]
	#print dist
	if(dist > 8):
		ax.annotate(row.name, list(row.values),
			xytext=(10,-5), 
			textcoords='offset points',
			size=10, 
			color='darkslategrey')

#tsne = sklearn.manifold.TSNE(n_components=2, learning_rate=50, n_iter=10000, n_iter_without_progress=1000)

pca = decomposition.PCA(n_components=40)
pca.fit(expression_table)

pca_expression_table = pd.DataFrame(pca.transform(expression_table), index=expression_table.index)

for learning_rate in [10,50,100,200,400,1000]:
	for perplexity in [5,8,10,20,30,50]:
		print "fitting tSNE, step:",step
		tsne = sklearn.manifold.TSNE(
			n_components = 2,
			learning_rate = learning_rate,
			perplexity = perplexity,
			n_iter = 10000,
			n_iter_without_progress = 1000
		)
		tsne_transformed_expression_2d = pd.DataFrame(tsne.fit_transform(pca_expression_table.values), index=pca_expression_table.index, columns=["x","y"])
		fig,ax = plt.subplots(figsize=(15, 10))
		ax = tsne_transformed_expression_2d.plot.scatter(x="x", y="y", s = annotations["day_size"].values, c=annotations["treatment_color"].values, ax=ax)
		ab=tsne_transformed_expression_2d.apply(annotate_df, axis=1)
		fn = "PCA_tSNE_training_lr_per/tSNE_lr."+str(learning_rate)+"_perplexity."+str(perplexity)+".png"
		plt.savefig(fn)
		plt.close()

tsne3d = sklearn.manifold.TSNE(n_components=3, learning_rate=20, n_iter=10000, n_iter_without_progress=1000, perplexity = 10)
tsne_transformed_expression_3d = pd.DataFrame(tsne3d.fit_transform(pca_expression_table.values), index=pca_expression_table.index, columns=["x","y","z"])
# plotly 3d plot
def plot_using_plotly(transformed_expression):
	import plotly.plotly as py
	import plotly.graph_objs as go
	layout = dict(
	width=1600,
	height=1080,
	autosize=False,
	#title='Test',
	scene=dict(
		xaxis=dict(
			gridcolor='rgb(0, 0, 0)',
			zerolinecolor='rgb(255, 0, 0)',
			showbackground=True,
			backgroundcolor='#bababa'
		),
		yaxis=dict(
			gridcolor='rgb(0, 0, 0)',
			zerolinecolor='rgb(255, 0, 0)',
			showbackground=True,
			backgroundcolor='#bababa'
		),
		zaxis=dict(
			#title="PC "+str(pc[2]+1),
			gridcolor='rgb(0, 0, 0)',
			zerolinecolor='rgb(255, 0, 0)',
			showbackground=True,
			backgroundcolor='#bababa'
		),
		#aspectratio = dict( x=1, y=1, z=0.7 ),
		aspectmode = 'manual'        
		),
	)
	data = []
	trace = dict(
		text = transformed_expression.index, #+ "\n" + transformed_expression["branch"],
		x = transformed_expression["x"],
		y = transformed_expression["y"],
		z = transformed_expression["z"],
		type = "scatter3d",    
		mode = 'markers',
		opacity = 0.80,
		marker = dict(
			size = annotations["day_size"].values,
			color=annotations["treatment_color"].values,
			symbol = annotations["symbol"]
			line=dict(width=1) )
		)
	data.append( trace )
	fig = dict(data=data, layout=layout)
	url = py.plot(fig, filename='single_cell-3d-tSNE', validate=False)

plot_using_plotly(tsne_transformed_expression_3d)
