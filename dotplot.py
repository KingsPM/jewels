#!/usr/bin/env python

__doc__ = '''
%prog query subject --qbed query.bed --sbed subject.bed
dotplot with lastz from 2 sequences
'''
__author__ = 'UNKOWN'


import os.path as op
import itertools
import sys
import collections
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle

sys.path.insert(0, op.join(op.dirname(__file__), ".."))
from grouper import Grouper

from dcbio.run.AlnRun import Align

class BlastLine(object):
    __slots__ = ('query', 'subject', 'pctid', 'hitlen', 'nmismatch', 'ngaps', \
                 'qstart', 'qstop', 'sstart', 'sstop', 'evalue', 'score', \
                 'qseqid', 'sseqid', 'qi', 'si')

    def __init__(self, sline):
        args = sline.split("\t")
        self.query = args[0]
        self.subject = args[1]
        self.pctid = float(args[2])
        self.hitlen = int(args[3])
        self.nmismatch = int(args[4])
        self.ngaps = int(args[5])
        self.qstart = int(args[6])
        self.qstop = int(args[7])
        self.sstart = int(args[8])
        self.sstop = int(args[9])
        self.evalue = float(args[10])
        self.score = float(args[11])

    def __repr__(self):
        return "BlastLine('%s' to '%s', eval=%.3f, score=%.1f)" % \
                (self.query, self.subject, self.evalue, self.score)

    def __str__(self):
        return "\t".join(map(str, [getattr(self, attr) \
                for attr in BlastLine.__slots__[:-4]]))

## BED parser
class Bed(list):

    def __init__(self, filename, key=None):
        self.filename = filename
        # the sorting key provides some flexibility in ordering the features
        # for example, user might not like the lexico-order of seqid
        self.key = key or (lambda x: (x.seqid, x.start, x.accn))
        for line in open(filename):
            if line[0] == "#": continue
            if line.startswith('track'): continue
            self.append(BedLine(line))

        self.seqids = sorted(set(b.seqid for b in self))
        self.sort(key=self.key)

    def get_order(self):
        return dict((f.accn, (i, f)) for (i, f) in enumerate(self))

    def get_simple_bed(self):
        return [(b.seqid, i) for (i, b) in enumerate(self)]

class BedLine(object):
    # the Bed format supports more columns. we only need
    # the first 4, but keep the information in 'stuff'.
    __slots__ = ("seqid", "start", "end", "accn", "stuff")

    def __init__(self, sline):
        args = sline.strip().split("\t")
        self.seqid = args[0]
        self.start = int(args[1])
        self.end = int(args[2])
        self.accn = args[3]
        self.stuff = args[4:] if len(args) > 4 else None

    def __str__(self):
        s = "\t".join(map(str, [getattr(self, attr) \
                    for attr in BedLine.__slots__[:-1]]))
        if self.stuff:
            s += "\t" + "\t".join(self.stuff)
        return s

    def __getitem__(self, key):
        return getattr(self, key)


## plotting functions
def get_breaks(bed):
    # get chromosome break positions
    simple_bed = bed.get_simple_bed()
    for seqid, ranks in itertools.groupby(simple_bed, key=lambda x:x[0]):
        ranks = list(ranks)
        # chromosome, extent of the chromosome
        yield seqid, ranks[0][1], ranks[-1][1]


def draw_box(clusters, ax, color="b"):

    for cluster in clusters:
        xrect, yrect = zip(*cluster)
        xmin, xmax, ymin, ymax = min(xrect), max(xrect), \
                                min(yrect), max(yrect)
        ax.add_patch(Rectangle((xmin, ymin), xmax-xmin, ymax-ymin,\
                                ec=color, fc='y', alpha=.5))
        #ax.plot(xrect, yrect, 'r.', ms=3)


def dotplot(segments, qbed, sbed, image_name, szoom):

    blast_fh = file(blast_file)
    blasts = [BlastLine(line) for line in blast_fh]
    seen = set()

    qorder = qbed.get_order()
    sorder = sbed.get_order()

    data = []
    for b in blasts:
        query, subject = b.query, b.subject
        if query not in qorder or subject not in sorder:
            continue
        key = query, subject
        if key in seen:
            continue
        seen.add(key)

        qi, q = qorder[query]
        si, s = sorder[subject]
        data.append((qi, si))

    fig = plt.figure(1,(8,8))
    root = fig.add_axes([0,0,1,1]) # the whole canvas
    ax = fig.add_axes([.1,.1,.8,.8]) # the dot plot

    x, y = zip(*data)
    ax.scatter(x, y, c='k', s=.05, lw=0, alpha=.9)

    xlim = (0, len(qbed))
    ylim = (0, len(sbed))

    xchr_labels, ychr_labels = [], []
    ignore = True # tag to mark whether to plot chromosome name (skip small ones)
    ignore_size = 100
    # plot the chromosome breaks
    for (seqid, beg, end) in get_breaks(qbed):
        ignore = abs(end-beg) < ignore_size
        xchr_labels.append((seqid, (beg + end)/2, ignore))
        ax.plot([beg, beg], ylim, "g-", alpha=.5)

    for (seqid, beg, end) in get_breaks(sbed):
        ignore = abs(end-beg) < ignore_size
        ychr_labels.append((seqid, (beg + end)/2, ignore))
        ax.plot(xlim, [beg, beg], "g-", alpha=.5)

    # plot the chromosome labels
    for label, pos, ignore in xchr_labels:
        pos = .1 + pos * .8/xlim[1]
        if not ignore:
            root.text(pos, .91, r"%s" % label, color="b", size=9, alpha=.5, rotation=45)

    for label, pos, ignore in ychr_labels:
        pos = .1 + pos * .8/ylim[1]
        if not ignore:
            root.text(.91, pos, r"%s" % label, color="b", size=9, alpha=.5, ha="left", va="center")

    # create a diagonal to separate mirror image for self comparison
    if is_self:
        ax.plot(xlim, ylim, 'm-', alpha=.5, lw=2)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # i always like the latex font
    _ = lambda x: r"$\mathsf{%s}$" % x.replace("_", " ").replace(" ", r"\ ")
    to_ax_label = lambda fname: _(op.basename(fname).split(".")[0])

    # add genome names
    ax.set_xlabel(to_ax_label(qbed.filename))
    ax.set_ylabel(to_ax_label(sbed.filename))

    # beautify the numeric axis
    [tick.set_visible(False) for tick in ax.get_xticklines() + ax.get_yticklines()]
    formatter = ticker.FuncFormatter(lambda x, pos: r"$%d$" % x)
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)
    plt.setp(ax.get_xticklabels() + ax.get_yticklabels(), color='gray', size=10)

    root.set_axis_off()
    print >> sys.stderr, "print image to %s" % image_name
    plt.savefig(image_name, dpi=1000)






if __name__ == "__main__":

    from Bio import SeqIO
    from itools import listpairs

    # readsequences
    records = []
    for seq_record in SeqIO.parse(sys.stdin, "fasta"):
        records.append(seq_record)

    # align
    for pair in listpairs(records):
        print >> sys.stderr, '\t', pair[0].id, pair[1].id
        alnrun = AlnRun.Align(pair)
        single = alnrun.lastz()

    # plot



'''

    import optparse

    parser = optparse.OptionParser(__doc__)
    parser.add_option("--qbed", dest="qbed", help="path to qbed")
    parser.add_option("--sbed", dest="sbed", help="path to sbed")
    parser.add_option("--format", dest="format", default="png", help="generate image of format (PNG, pdf, ps, eps, svg, etc.))
    (options, args) = parser.parse_args()

    if not (len(args) == 2 and options.qbed and options.sbed):
        sys.exit(parser.print_help())

    qbed = Bed(options.qbed)
    sbed = Bed(options.sbed)

    blast_file = args[0]

    image_name = op.splitext(blast_file)[0] + "." + options.format
    dotplot(blast_file, qbed, sbed, image_name)

'''