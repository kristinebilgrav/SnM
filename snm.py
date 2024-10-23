import pysam
import sys

def import_haploblocks(file):
    dict={}    
    for line in open(file):
        chr=line.split('\t')[0]
        start=int(line.split('\t')[1])
        end=int(line.split('\t')[2])
        id=line.rstrip('\n').split('\t')[3]
        if chr not in dict:
            dict[chr]={}
        if id not in dict[chr]:
            dict[chr][id]=(start, end)

    return dict


def coverage(block_dict, bam):
    samfile = pysam.AlignmentFile(bam, "rb" )
    outlist=[]

    for chr in block_dict:
        for blockID in block_dict[chr]:
            start= block_dict[chr][blockID][0]
            end= block_dict[chr][blockID][1]
            total=[]
            covhap1=0
            covhap2=0
            covnohap=0
            for pileupcolumn in samfile.pileup(chr, start, end): #pileupcolumn: all the reads that map to a single base in the region
                #print(pileupcolumn)
                
                total.append(pileupcolumn.n)
                #print('total', total)
                for pileupread in pileupcolumn.pileups:
                    try:
                        haplotype = pileupread.alignment.get_tag('HP')
                    except:
                        haplotype = 0
                    if haplotype == 1:
                        covhap1 += 1
                    elif haplotype == 2:
                        covhap2 +=1
                    else:
                        covnohap +=1


            ratio=covhap1/covhap2 # Output ratio to see copynumber
            if covhap1 > covhap2:
                dominant= 'hap1'
            elif covhap2 > covhap1:
                dominant='hap2'
            else:
                dominant='unphased'

            avrgcov=sum(total)/len(total)
            lst= [chr, str(start), str(end), blockID, str(covhap1), str(covhap2), str(covnohap), str(ratio), dominant,str(avrgcov)]    

            outlist.append(lst)

    return outlist



blocks= import_haploblocks(sys.argv[1])

res=coverage(blocks, sys.argv[2])

output = open(sys.argv[3], 'w')
for l in res:
    output.write('\t'.join(l) + '\n')
