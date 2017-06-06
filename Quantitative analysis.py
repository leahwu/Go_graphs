import plotly.plotly as py
from plotly.tools import FigureFactory as FF
from DCMGenerator import *

n = 2000

spear_corr_matrix = [["Correlation"]+["d = "+str(1 + 0.2*i) for i in range(5)]]
for j in range(16):
    spear_corr_matrix.append(["beta = "+ str(2.5 + 0.5*j)]+ [DCMGenerator(1, 1+0.2 * i, 2.5+0.5*j, n,'Erased').spearman_test()[0] for i in range(5)])

corr_table = FF.create_table(spear_corr_matrix, index=True)
for i in range(len(corr_table.layout.annotations)):
    corr_table.layout.annotations[i].font.size = 20

py.iplot(corr_table, filename='corr_table')