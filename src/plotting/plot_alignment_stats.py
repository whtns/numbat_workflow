#!/usr/bin/python 

# import system library - takes care of reading arguments from command line
import sys
# pandas library helps with reading csv file and with statistics
import pandas as pd
# library matplotlib takes care of plotting the graphs
import matplotlib.pyplot as plt
# library saborn provides more pleasant color palletes for the graphs
import seaborn as sns
import re
# argv[1] = tab delimited csv file with gene counts and percentages such as:
#	Vanilla from sequencer	After trimming	IN:reads	IN:paired	%	PAIR-AL-C:0	%	PAIR-AL-C:1	%	PAIR-AL-C:>1	%	PAIR-AL-D:1	%	MATE-AL:0	%	MATE-AL:1	%	MATE-AL:>1	%	IN:unpaired	%	AL:0	%	AL:1	%	AL:>1	%	overall alignment rate
#1	9026440	8603753	4507192	4096561	90.89	137198	3.35	3548995	86.63	410368	10.02	21287	15.52	124769	53.82	57514	24.81	49539	21.37	410631	9.11	14982	3.65	316268	77.02	79381	19.33	98.38
#2	10189312	9754053	5087553	4666500	91.72	141030	3.02	4019108	86.13	506362	10.85	20894	14.82	122301	50.9	62563	26.04	55408	23.06	421053	8.28	14820	3.52	321013	76.24	85220	20.24	98.59

# argv[2] = file with list of cells to plot
# 100
# 101
# 102



sheet_file = sys.argv[1]
cells_to_display_file = sys.argv[2]
output_file = re.sub(".csv", "", sheet_file)
cells = open(cells_to_display_file).read().splitlines()
bar_figsize = (10,6)
dpi=400

columns_with_pairs = ["Vanilla from fastqc","IN:reads","IN:paired", "PAIR-AL-C:0", "PAIR-AL-C:1", "PAIR-AL-C:>1", "PAIR-AL-D:1"]

data = pd.read_csv(sheet_file, sep="\t", index_col=0)
data.index = data.index.astype(str)

data = data.loc[cells,]

data[columns_with_pairs] = data[columns_with_pairs]*2
print data.head()
data = data.sort_values(by="Vanilla from fastqc")

columns = []
columns.append(["Vanilla from fastqc"])
columns.append(["IN:reads"])
columns.append(["PAIR-AL-C:1", "PAIR-AL-C:>1", "PAIR-AL-D:1", "MATE-AL:1", "MATE-AL:>1", "AL:1", "AL:>1"]) #  "AL:1", "AL:>1"


pos =    1.0
offset = 1.0
width  = 0.45
fig, ax = plt.subplots(1,figsize=bar_figsize)
ax = data[columns[0]].plot.bar(stacked=True, position=pos, width=width, ax=ax, color="blue")
ax = data[columns[1]].plot.bar(stacked=True, position=pos, width=width, ax=ax, color="green")
ax = data[columns[2]].plot.bar(stacked=True, position=pos-offset, width=width, ax=ax, colormap="Set2")
plt.tick_params(axis='both', which='major', labelsize=1)
plt.title(sheet_file)

x1,x2 = ax.get_xlim()
ax.set_xlim(-0.5,x2)
fig.tight_layout()
#~ plt.show()


fig.savefig(output_file+".pdf", dpi=dpi)


