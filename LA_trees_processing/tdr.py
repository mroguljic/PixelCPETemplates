from scipy import stats
import math
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import statsmodels.tools.numdiff as nd
import sys
# argument is the name of the file to open
if len(sys.argv) != 2:
  sys.exit(f'wrong number of arguments = {len(sys.argv)-1}')
else:
# open the data file for reading only
  f = open(sys.argv[1])
#define lists to contain the data sets and data
xset = []
yset = []
yeset = []
cset = []
lset = []
lsset = []
fmtset = []
xdata = []
ydata = []
yedata = []
def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False
# loop through the file and parse the data elements on each line
nlines = 0
for line in f:
# skip blank lines and those that begin with #
     words = line.split()
     if len(words) < 1 or line[:1] == "#":
        continue
# the first two are the x- and y-axis labels
     nlines += 1
     if nlines == 1:
        xtitle = line
     elif nlines == 2:
        ytitle = line
     else:
        if isfloat(words[0]):
           x = float(words[0])
           y = float(words[1])
           xdata.append(x)
           ydata.append(y)
           if len(words) > 2:
              ye = float(words[2])
              yedata.append(ye)
        else:
           lset.append(words[0])
           cset.append(words[1])
           xd = np.array(xdata)
           yd = np.array(ydata)
           yed = np.array(yedata)
           yset.append(yd)
           xset.append(xd)
           yeset.append(yed)
           xdata.clear()
           ydata.clear()
           yedata.clear()
           if len(words) < 3:
              lsset.append("solid")
              fmtset.append("o")
           else:
              lsset.append(words[2])
              if len(words) < 4:
                 fmtset.append("o")
              else:
                 fmtset.append(words[3])
           
# don't forget to close the file
f.close

nsets = len(lset)

for i in range(nsets):
   npoints = len(xset[i]) # we'll need the number of points to plot the fit function
   print(f'data set {i} with {npoints} points')
plt.figure(constrained_layout=True)
plt.title(sys.argv[1], fontsize=16)
plt.xlabel(xtitle,fontsize=16)
plt.ylabel(ytitle,fontsize=16,labelpad=-10)
# plt.tick_params(bottom="on", left="on")
for i in range(nsets):
# plot each set
   if len(yeset[i]) == 0:
     plt.plot(xset[i], yset[i], color=cset[i], linestyle=lsset[i], label=lset[i])
   else:
     plt.errorbar(xset[i], yset[i], yerr=yeset[i], fmt=fmtset[i], linestyle=lsset[i], color=cset[i], label=lset[i])
# label the plot
plt.ylim(bottom=0)
plt.legend()
outfile = 'data_mc_plots/'+sys.argv[1] + '.pdf'
plt.savefig(outfile)
plt.show()
