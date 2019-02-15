#!/usr/bin/env python


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

    def disable(self):
        self.BLACK = ''
        self.CYAN = ''
        self.WHITE = ''
        self.PURPLE = ''
        self.BLUE = ''
        self.GREEN = ''
        self.YELLOW = ''
        self.RED = ''
