# -*- coding: UTF-8 -*-
#!/opt/local/bin/python

import os
import sys
import json
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pandas as pd
#-----------------------------
scriptpath = os.path.dirname(os.path.realpath(__file__))
defaultconfig = os.path.join(scriptpath, "configs/gwas.json")
#-----------------------------
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--scorefile',  default="/Users/jkb/Desktop/temp/metal-semode-QCspecific-filteredd05AF25.tsv.sorted.tsv",  help='p value file')
# /mnt/ris-fas1a/linux_users2/kbaillie/crispr/screens_corrected/rawdata/combineds1s2.correct.txtcombinedscreens_corrected_P6_vs_P4.qnz
parser.add_argument('-fi', '--scorefile_id_col', default='gene',    help='column name')
parser.add_argument('-fs', '--scorefile_score_col', default='P',    help='column name - score to be plotted')
parser.add_argument('-ft', '--scorefile_score_type', default='p', choices=['p','z','logp', 'fdr'])
parser.add_argument('-sfs', '--scorefilesep', default='\t')
parser.add_argument('-pfs', '--posfilesep', default='\t')
parser.add_argument('-t', '--chart_title', default='chart')
parser.add_argument('-s', '--siglines',    action='append', default=[], help='add as many significance threshold lines as you like')
parser.add_argument('-ha', '--highlight_all',    action="store_true", default=False,    help='plot red spots for all guides in highlightlist')
parser.add_argument('-nc', '--numcolors', type=int, default=2)
parser.add_argument('-clf', '--chromlengthsfile', default='null',    help='specify chrom lengths')
parser.add_argument('-p', '--positionfile', default="use_scorefile",    help='a file that lists the positions on genome of each identifier')
parser.add_argument('-pi', '--positionfile_id_col', default='gene',    help='column name')
parser.add_argument('-pc', '--positionfile_chrom_col', default='CHR',    help='column name')
parser.add_argument('-pa', '--positionfile_address_col', default='BP',    help='column name')
parser.add_argument('-hc', '--highlight_column', default="",    help='column name to choose hits to highlight')
parser.add_argument('-ht', '--highlight_threshold', default="",    help='threshold in highlight_column')
parser.add_argument('-htt', '--highlight_test', default="less",    help='less or more. less would mean that everything below highlight_threshold is red, more would mean above')
parser.add_argument('-c', '--configfile', default=defaultconfig,    help='config file')
parser.add_argument('-o', '--outputformat', default='png',    help='pdf png or svg')
parser.add_argument('-ss', '--skipsex',    action="store_true", default=False,    help='skip sex and MT etc')
parser.add_argument('-ms', '--scorefile_miami_col', default='null',    help='column name - score to be plotted on negative side')
parser.add_argument('-cs', '--chrom_spacing', default=25000000,  type=int,  help='column name - score to be plotted on negative side')
args = parser.parse_args()
#-----------------------------
boxcolor = "#f4f4f4"
point_size = 2
sig_point_size = 2
sig_color = 'red'
axis_text_size = 10
chr_text_size = 8
dpi = 300
xlab = "Chromosome"
ylab = "Z score"
if args.scorefile_score_type == 'p' or args.scorefile_score_type == 'logp':
    ylab = r'$-\log_{10}$ (p value)'
if args.scorefile_score_type == 'fdr':
    ylab = r'$-\log_{10}$ (FDR)'
#-----------------------------
def infiniteindex(i, thislen):
    '''
        correct a list index that is greater than the length of the list
    '''
    thislen = float(thislen)
    i = float(i)
    x = int( round(thislen*(i/thislen - int (i/thislen)),4  ))
    return x

from palettable.colorbrewer.qualitative import Set1_9
from palettable.colorbrewer.qualitative import Set2_8
from palettable.colorbrewer.qualitative import Set3_12

def getcolor(thisinteger, maxcolors=12):
    try:
        configdata['plotcolors']
    except:
        return Set2_8.mpl_colors[infiniteindex(thisinteger, maxcolors)]
    return configdata['plotcolors'][infiniteindex(thisinteger, len(configdata['plotcolors']))]

#-----------------------------
with open(args.configfile) as json_data_file:
    configdata = json.load(json_data_file)
graphout = os.path.expanduser(configdata["graphdir"])
datasetlabel = configdata["label"]
highlightlist = {}
for key in configdata["manhattan_highlightlist"]:
    '''
    input can be:
        "GENE" or
        ["GENE",2,2] or
        "identifier":["GENE",2,2] or
        "identifier":"GENE"
    '''
    try:
        value = configdata["manhattan_highlightlist"][key]
        # now we know that input is "identifier":["GENE",2,2] or "identifier":"GENE"
    except:
        # now input could be either "GENE" or ["GENE",2,2]
        value = key
        if type(key) is list:
            key = key[0]
    # now we know that key is "identifier" and value is "GENE" or ["GENE",2,2]
    if type(value) is not list:
        value = [value, 1,1]
    highlightlist[key] = value

#-----------------------------
chromswapper = [
    ["X",23],
    ["Y",24],
    ["XY",25],
    ["MT",26],
    ["M",26],
]
chrletter2num = {x[0]:x[1] for x in chromswapper}
chrnum2letter = {x[1]:x[0] for x in chromswapper}
#-----------------------------

sdf = pd.read_table(args.scorefile, sep=args.scorefilesep) # scores data frame
if args.positionfile == "use_scorefile":
    df = sdf
else:
    pdf = pd.read_table(args.positionfile, sep=args.posfilesep) # positions data frame
    print (pdf)
    df = pd.merge(sdf, pdf, how='inner', left_on=args.scorefile_id_col, right_on=args.positionfile_id_col)
print (df)
print (df.shape)
if len(args.highlight_column)>0 and len(args.highlight_threshold)>0:
    df["highlight"] = df[args.highlight_column]
else:
    df["highlight"] = False
if args.scorefile_miami_col == "null":
    df = df[[args.scorefile_id_col, args.positionfile_chrom_col, args.positionfile_address_col, args.scorefile_score_col, "highlight"]]
    df.columns = ['gene', 'chrom', 'position', 'score', 'highlight']
else:
    df = df[[args.scorefile_id_col, args.positionfile_chrom_col, args.positionfile_address_col, args.scorefile_score_col, args.scorefile_miami_col,  "highlight"]]
    df.columns = ['gene', 'chrom', 'position', 'score', 'miami_score', 'highlight']

df.loc[:,'chrom'] = df["chrom"].astype("str").str.split("_", n = 1, expand = True)[0].str.replace('chr','')
df.loc[:,'chrom'] = df['chrom'].replace(chrletter2num)
df.loc[:,'chrom'] = pd.to_numeric(df['chrom'], errors='coerce')
df = df.dropna() # drop all rows with any nans
df.drop_duplicates(keep='first',inplace=True)

if args.scorefile_score_type == 'p' or args.scorefile_score_type == 'fdr':
    df['score'] = -np.log10(df['score'])
    if args.scorefile_miami_col != "null":
        df['miami_score'] = -np.log10(df['miami_score'])
        print (df[['score', 'miami_score']])

customchromlengths = {}
if args.chromlengthsfile != 'null':
    with open(args.chromlengthsfile) as f:
        lines = [x.strip().split(": ") for x in f.readlines() if len(x)>1]
        for x in lines:
            x[0] = x[0].replace('chr','')
            try:
                c = int(chrletter2num[x[0]])
            except:
                c = int(x[0])
            customchromlengths[c] = int(x[1])

# make plot
plt.ioff()
figure = None
figure = plt.figure(figsize=(18, 6), frameon=True)
#figure = plt.figure(figsize=(12, 6), frameon=True)

# The chromosome spacing
centimorgans = False
if centimorgans:
    args.chrom_spacing = 25.0

# Creating the ax and modify it
ax = figure.add_subplot(111)
ax.xaxis.set_ticks_position("none")
ax.yaxis.set_ticks_position("left")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["bottom"].set_visible(False)

ax.set_xticks([])
ax.set_xticklabels([])
ax.set_ylabel(ylab)
ax.set_xlabel(xlab)
ax.set_title(args.chart_title, fontsize=16)

# Now plotting for each of the chromosome
starting_pos = 0
annots = []
ticks = []

allchroms = list(df['chrom'].unique()) + list(customchromlengths.keys())
print (allchroms)
allchroms = sorted(list(set(allchroms)))
found = []
for i, chrom in enumerate(allchroms):
    if args.skipsex and chrom > 22:
        continue

    color = getcolor(i, args.numcolors)
    cd = df.loc[df['chrom'] == chrom] # chrom data

    try:
        max_pos = customchromlengths[chrom]
    except:
        max_pos = cd['position'].max()

    # The box
    xmin = starting_pos - (args.chrom_spacing / 2)
    xmax = max_pos + starting_pos + (args.chrom_spacing / 2)
    if i % 2 == 1:
        ax.axvspan(xmin=xmin, xmax=xmax, color=boxcolor)

    # The chromosome label
    ticks.append((xmin + xmax) / 2)

    ax.plot(cd['position'] + starting_pos, cd['score'],
            marker="o", ms=point_size, mfc=color, mec=color,
            ls="None")

    if args.scorefile_miami_col != "null":
        ax.plot(cd['position'] + starting_pos, -cd['miami_score'],
                marker="o", ms=point_size, mfc=color, mec=color,
                ls="None")

    for identifier in highlightlist:
        if identifier in cd['gene'].values:
            found.append(identifier)
            print ("identifier found", identifier)
            thisgene = cd[cd['gene']==identifier].reset_index()
            if args.highlight_all:
                ax.plot(thisgene['position'] + starting_pos, thisgene['score'],
                marker="o", ms=sig_point_size, mfc=sig_color, mec=sig_color,
                ls="None")
            if len(thisgene)>1:
                thisgene = thisgene.iloc[thisgene['score'].idxmax()] # needed because some genes have multiple hits. take the best score.
            spotcoords = [float(x) for x in [thisgene['position'] + starting_pos, float(thisgene['score'])]]
            x_displacement = 30000000*highlightlist[identifier][1]
            y_displacement = 0.5*highlightlist[identifier][2]
            print (identifier, highlightlist[identifier], x_displacement, y_displacement)
            if float(thisgene['score'])>0:
                textcoords = [spotcoords[0]+x_displacement, spotcoords[1]+y_displacement]
                l = mlines.Line2D([spotcoords[0],textcoords[0]-50000], [spotcoords[1], textcoords[1]],  linestyle='-', linewidth=1.2, color='#3a3a3a') #color='#3a3a3a'
            else:
                textcoords = [spotcoords[0]+x_displacement, spotcoords[1]-y_displacement]
                l = mlines.Line2D([spotcoords[0],textcoords[0]-50000], [spotcoords[1], textcoords[1]+0.3],  linestyle='-', linewidth=1.2, color='#3a3a3a') #color='#3a3a3a'
            font = {
                'family': 'serif',
                'color':  'black',
                #'weight': 'bold',
                'style' : 'italic'
                }
            ax.text(textcoords[0], textcoords[1], highlightlist[identifier][0], fontdict=font)
            ax.add_line(l)

    if len(args.highlight_column)>0 and len(args.highlight_threshold)>0:
        if args.highlight_test == "less":
            hcd = cd[cd["highlight"] < float(args.highlight_threshold)]
        elif args.highlight_test == "more":
            hcd = cd[cd["highlight"] > float(args.highlight_threshold)]
        else:
            print ("{} not understood".format(args.highlight_test))
        ax.plot(hcd['position'] + starting_pos, hcd['score'],
                marker="o", ms=point_size, mfc=sig_color, mec=sig_color,
                ls="None")

    # Changing the starting point for the next chromosome
    starting_pos = max_pos + starting_pos + args.chrom_spacing

found = set(found)
sought = set(highlightlist.keys())
missed = sought - found
if len(missed)>0:
    print ("\n**You looked for these identifiers but they weren't found:**")
    print ("\n".join(list(missed)))
    print ("\n")

ax.set_xlim(0 - args.chrom_spacing, starting_pos + args.chrom_spacing)
if len(configdata['ylim'])>0:
    ax.set_ylim(configdata['ylim'][0], configdata['ylim'][1])

# Putting the xticklabels
ax.set_xticks(ticks)
newlabels = []
for c in allchroms:
    try:
        c = int(c)
    except:
        pass
    try:
        newlabels.append(chrnum2letter[c])
    except:
        newlabels.append(c)
ax.set_xticklabels(newlabels)

for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(axis_text_size)

for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(chr_text_size)

for abline_position in args.siglines:
    a = float(abline_position)
    ax.axhline(a, color="red", ls="--", lw=1.2)

# Saving or plotting the figure
mpl.rcParams['savefig.dpi'] = dpi
mpl.rcParams['ps.papersize'] = "auto"
mpl.rcParams['savefig.orientation'] = "landscape"
outfilename = '{}.{}'.format(".".join([datasetlabel, args.scorefile_score_type, args.chart_title]),args.outputformat).replace(" ","_")
figfile = os.path.join(graphout, outfilename)
print (figfile)
plt.savefig(figfile, bbox_inches="tight")




