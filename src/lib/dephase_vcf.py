
from tempfile import mkstemp
import shutil
import gzip
import re
from os import fdopen, remove
import os

def dephase_vcf(vcf_in, loc):
    """
    Adapted from @apoursh
    """
    def replacemany(adict, astring):
        pat = '|'.join(re.escape(s) for s in adict)
        there = re.compile(pat)
        def onerepl(mo): return adict[mo.group()]
        return there.sub(onerepl, astring)

    ext = vcf_in.split(os.extsep)
    outname = os.path.join(loc, ext[0] + '.dephased.' + '.'.join(ext[1:]))
    if os.path.isfile(outname): #check if it already exists
        return
    rep = {"0|1": "1/0", "0|0": "0/0", "1|0":"1/0", "1|1":"1/1"}
    #Produce temp file
    fh, abs_path = mkstemp()
    with fdopen(fh,'w') as new_file:
        with gzip.open(vcf_in, 'rt') as old_file:
            for line in old_file:
                if line[0] != '#':
                    line = replacemany(rep, line)
                new_file.write(line)
    #Move temp file to new location
    with open(abs_path, 'rb') as f_in:
        with gzip.open(outname, 'wb') as loc:
            print ("Copying dephased vcf over to: {}".format(loc))
            shutil.copyfileobj(f_in, loc)
    remove(abs_path)

    return outname