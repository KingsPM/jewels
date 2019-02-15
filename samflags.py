#!/usr/bin/env python


import sys

class Colors:
    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    PURPLE = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'

    BGBLACK = '\033[40m'
    BGRED = '\033[41m'
    BGGREEN = '\033[42m'
    BGYELLOW = '\033[43m'
    BGBLUE = '\033[44m'
    BGPURPLE = '\033[45m'
    BGCYAN = '\033[46m'
    BGWHITE = '\033[47m'

    RESET = '\033[0m'
    OFF = '\033[0m'
    BOLD = '\033[1m'
    UNDERSCORE = '\033[4m'
    BLINK = '\033[5m'
    REVERSE = '\033[7m'
    CONCEAL = '\033[8m'

def help():
    print >> sys.stderr, '''
    samflags.py +[requires flag] -[absent flag] < <SAMFILE>

    FLAG BIT    DESCRIPTION
    ---- -----  ------------------------------------------------------
       1 0x1    template having multiple segments in sequencing
       2 0x2    each segment properly aligned according to the aligner
       4 0x4    segment unmapped
       8 0x8    next segment in the template unmapped
      16 0x10   SEQ being reverse complemented
      32 0x20   SEQ of the next segment in the template being reversed
      64 0x40   the Frst segment in the template
     128 0x80   the last segment in the template
     256 0x100  secondary alignment
     512 0x200  not passing quality controls
    1024 0x400  PCR or optical duplicate
    '''

if __name__=="__main__":
    # only return when flag set (for generating ad-hoc mapping statistics)
    try:
        askFlags = map(int,sys.argv[1:])
    except:
        help()
        sys.exit(1)

    requiredFlags, notFlags = [], []
    for fl in askFlags:
        if fl < 0:
            notFlags.append(abs(fl))
        else:
            requiredFlags.append(fl)

    for line in sys.stdin:
        if line.startswith('@'):
            continue
        f = line.split()

        # check if passes flag requirements
        if (requiredFlags and not all([ int(f[1]) & x for x in requiredFlags ])) or (notFlags and any([ int(f[1]) & x for x in notFlags ])):
            continue

        # build flag array
        flags = []
        if int(f[1]) & 1:
            flags.append('2') # MULTIPLE SEGMENTS
        else:
            flags.append('1')
        if int(f[1]) & 2:
            flags.append(Colors.GREEN+'*'+Colors.RESET) # PROPER PAIR
        else:
            flags.append(Colors.RED+'x'+Colors.RESET)


        def active():
            if not int(f[1]) & 4:  # mapped
                if int(f[1]) & 16:
                    return Colors.BOLD+Colors.GREEN+'<'+Colors.RESET # REVERSE
                else:
                    return Colors.BOLD+Colors.GREEN+'>'+Colors.RESET # REVERSE
            else:  # unmapped
                return Colors.RED+'-'+Colors.RESET

        def inactive():
            if not int(f[1]) & 8:  # mapped
                if int(f[1]) & 32:
                    return '<' # REVERSE
                else:
                    return '>'
            else:
                return '-'

        # POSITION AND ORIENTATION
        if int(f[1]) & 64: # First segment on the left or last
            if f[8].startswith('-'):
                flags.append(inactive())
                flags.append(active())
            else:
                flags.append(active())
                flags.append(inactive())

        elif int(f[1]) & 128:  # Last segment
            if f[8].startswith('-') or int(f[8]) == 0:
                flags.append(inactive())
                flags.append(active())
            else:
                flags.append(active())
                flags.append(inactive())

        # MAPPING
        if int(f[1]) & 256:
            flags.append(Colors.CYAN+'2'+Colors.RESET) # SECONDARY
        else:
            flags.append(Colors.BLUE+'1'+Colors.RESET)
        if int(f[1]) & 512:
            flags.append(Colors.RED+'B'+Colors.RESET) # BAD QUAL
        else:
            flags.append(Colors.GREEN+'G'+Colors.RESET) # GOOD QUAL
        if int(f[1]) & 1024:
            flags.append(Colors.RED+'D'+Colors.RESET) # DUPLICATE
        else:
            flags.append(' ')
        if int(f[1]) & 2048:  # supplementary alignment
            flags.append(Colors.YELLOW+'+'+Colors.RESET) # DUPLICATE
        else:
            flags.append(' ')

        print '\t'.join([''.join(flags)] + f[0:1] + f[2:9])

