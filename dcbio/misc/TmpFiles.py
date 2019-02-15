#!/usr/bin/env python

''' simple helper class for handling temporary files'''
import sys
import string

class TmpFiles:
    def __init__(self,number=0):
        self.handles = []
        for i in range(number):
            self._tmpfile()
        return

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        for fh in self.handles:
            try:
                fh.close()
            except:
                pass
            else:
                sys.stderr.write('(TmpFile: closed '+ fh.name + ')')
            # remove file
            try:
                os.unlink(fh.name)
            except:
                pass
            else:
                sys.stderr.write('(TmpFile: unlinked '+ fh.name + ')')
        return

    def _tmpfile(self):
        self.handles.append(tempfile.NamedTemporaryFile(delete=False))
        return self.handles[-1]


class Filename(str):
    def __init__(self,filename):
        self = filename 
        return

    def subst(self,s='_'):
        '''filename with substituted characters'''        
        valid_chars = "-_.()/%s%s" % (string.ascii_letters, string.digits)
        fn = ''.join([ x if x in valid_chars else s for x in self ])
        try:
            assert len(fn) > 0
        except:
            raise Exception('Filename has no valid characters')
        return fn


if __name__=="__main__":
    fi = Filename("This Is a (valid) - filename%$&$ .txt")
    print 'FI', fi
    print 'MH', fi.subst()
   
