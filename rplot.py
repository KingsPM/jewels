#!/usr/bin/env python

from rpy2.robjects.vectors import DataFrame
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import *
#import rpy2
#print dir(rpy2.robjects.vectors)

#mtcars = datasets.__rdata__.fetch('mtcars')['mtcars']

def scatter_plot(data,x,y,size,color):
    import math, datetime
    import rpy2.robjects.lib.ggplot2 as ggplot2
    import rpy2.robjects as ro
    base = importr('base')
    grdevices.pdf(outfile)
    gp = ggplot2.ggplot(data)
    pp = gp + ggplot2.aes_string(x=x, y=y, size=size, col=color, shape='factor(gear)') + ggplot2.geom_point()
    pp.plot()
    grdevices.dev_off()
    return

def plot_datalines(data_points, outfile, xlabels = None, title = ""):
    graphics = importr('graphics')
    grdevices = importr('grDevices')
    r_base = importr('base')
    grdevices.pdf(outfile)
    graphics.plot(data_points[0], data_points[1], type="o", xlab = '', ylab='', main=title)
    if xlabels != None:
        xlab = StrVector(xlabels)
        graphics.axis(1, at = IntVector(range(1,len(xlabels)+1)), lab=xlab)
        graphics.axis(2)
    grdevices.dev_off()
    return

def plot_barplot(data_points, outfile, xlabel = '', ylabel = '', title = "", line_width=2):
    graphics = importr('graphics')
    grdevices = importr('grDevices')
    r_base = importr('base')


    grdevices.postscript(outfile)

    X = FloatVector(data_points[0])
    Y = FloatVector(data_points[1])

    graphics.plot(X, Y, xlab=xlabel, ylab=ylabel, type='n', main = title)
    graphics.segments(X, Y, X,FloatVector([0]*len(data_points[0])), lwd=line_width)
    grdevices.dev_off()

    return


def plot_multiple_datalines(data_points, legend_col, outfile, xlabels=None, title="", xlabel='', ylabel=''):
    graphics = importr('graphics')
    grdevices = importr('grDevices')
    r_base = importr('base')
    grdevices.pdf(outfile)

   # graphics.par(new=True)
    for i in range(len(data_points)):
        data = data_points[i]
        graphics.plot(data[0], data[1], axes=False, type="o", col=legend_col[1][i], xlab=xlabel, ylab=ylabel, main=title)
        graphics.par(new=True)
    print legend_col[1]
    graphics.legend("topright", leg=StrVector(legend_col[0]), fill=StrVector(legend_col[1]), cex=0.5, bty="n")
    if xlabels != None:
        xlab = StrVector(xlabels)
        graphics.axis(1, at=IntVector(range(1, len(xlabels)+1)), lab=xlab)
        graphics.axis(2)
    grdevices.dev_off()

def plot_histogram(data_points, outfile, params, title="", xlabel='', ps=False, stats=True):
    # hist(mtcars$mpg, breaks=12, col="red")
    graphics = importr('graphics')
    grdevices = importr('grDevices')
    r_base = importr('base')

    minx, maxx, interval = params

    stats = importr('stats')
    if not ps:
        grdevices.pdf(outfile)
    else:
        grdevices.postscript(outfile)

    data = FloatVector(data_points)

    graphics.hist(data, r_base.seq(minx, maxx, interval), prob=True, col="gray", main=title, xlab=xlabel)
    if stats:
        graphics.lines(stats.density(data, bw=interval), col="red")

    graphics.rug(data)
    grdevices.dev_off()


if __name__ == "__main__":
    import random
    plot_histogram([random.randint(0, 5) for i in range(1000)], "testplot.pdf", (0, 5, 1), 'density of random integers')
