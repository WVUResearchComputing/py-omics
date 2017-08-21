#!/usr/bin/env python

import json

def process_gtf(inputpath, outputpath):

    rf=open(inputpath)
    ret=[]
    mRNA = None

    for line in rf.readlines():
        sline=line.split()
        if sline[2]=='mRNA':
            data={'chr': sline[0], 'gene': sline[2], 'start': int(sline[3]), 'stop': int(sline[4]), 'strand': sline[6], 'gene_id': sline[9][1:-2], 'transcript_id': sline[11][1:-2], 'regions':[]}
            if mRNA is not None:
                ret.append(mRNA)
                print(mRNA)
            mRNA=data
        else:
            data={'gene': sline[2], 'start': int(sline[3]), 'stop': int(sline[4])}
            mRNA['regions'].append(data)
    ret.append(mRNA)

    wf=open(outputpath, 'w')
    json.dump(ret, wf,sort_keys=True, indent=4, separators=(',', ': '))
    wf.close()


def transcrip_map(gtfpath, inputpath, outputpath):
    rf=open(gtfpath)
    data=json.load()

    rf=open(inputpath)
    rf.readline()

    ret={}
    for line in rf.readlines():

        sline=line.split()
        data2={'chr': sline[0], 'start': int(sline[1]), 'stop': int(sline[2]), 'perc_control': float(sline[5]), 'perc_estrogen': float(sline[6]), 'delta': float(sline[7])}

        #print('Searching for %s' % data2)
        for imRNA in data:
            if data2['chr'] == imRNA['chr']:
                toprint=False
                if data2['start'] >= imRNA['start'] and data2['stop'] <= imRNA['stop']:
                    inexon = False
                    overlap = False
                    for ireg in imRNA['regions']:
                        if data2['start'] >= ireg['start'] and data2['stop'] <= ireg['stop']:
                            inexon = True
                            overlap='Exon'
                            print('Start inside, Stop inside')
                        elif data2['start'] < ireg['start'] and data2['stop'] >= ireg['start']:
                            inexon=True
                            overlap='Intron-Exon'
                            print('Start before, Stop inside')
                        elif data2['start'] < ireg['stop'] and data2['stop'] >= ireg['stop']:
                            inexon=True
                            overlap='Exon-Intron'
                            print('Start inside, Stop after')
                    if inexon is False:
                        print("Intron: [%d - %d] is inside mRNA in [%d - %d]" % (data2['start'],data2['stop'], imRNA['start'], imRNA['stop']))
                    else:
                        print("%s: [%d - %d] is inside mRNA in [%d - %d]" % (overlap,data2['start'],data2['stop'], imRNA['start'], imRNA['stop']))
                    toprint=True
                elif data2['start'] < imRNA['start'] and data2['stop'] >= imRNA['start'] :
                    print("5'Overlap: [%d - %d] is inside mRNA in [%d - %d]" % (data2['start'],data2['stop'], imRNA['start'], imRNA['stop']))
                    toprint=True
                elif data2['start'] < imRNA['stop'] and data2['stop'] >= imRNA['stop'] :
                    print("3'Overlap: [%d - %d] is inside mRNA in [%d - %d]" % (data2['start'],data2['stop'], imRNA['start'], imRNA['stop']))
                    toprint=True
            if toprint:
                if imRNA['gene_id'] not in ret:
                    ret[imRNA['gene_id']]=[]
                ret[imRNA['gene_id']].append(data2['delta'])
                break

    wf=open('genes.json','w')
    json.dump(ret, wf,sort_keys=True, indent=4, separators=(',', ': '))
    wf.close()
