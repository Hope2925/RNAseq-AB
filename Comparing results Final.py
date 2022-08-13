
"""___________________________________________________________________________________________________________________
This script will go through the fasta & gff3 files of the A1S genome of AB and combined MFF genome and pAB3 genome of
AB and create a dictionary determining which A1S loci correspond to ACX60. It will then take expression files and
implement BLAST data (RefSeq protein ID, new locus, and BLAST results) and change the logFC to be more intuitive (-/+).

Date completed: December 3, 2022
Programmer Name: Hope Townsend (Kirby) (kirbyha)
____________________________________________________________________________________________________________________"""
import pandas as pd
import itertools




class Genome(object):
    def __init__(self, GFF3file, Fastafile, Version):
        """This function initializes appropriate dictionaries, values, and lists for the GFF3 file of study."""
        self.GFF_file = open(GFF3file, 'r')
        self.FASTA_file = open(Fastafile, 'r')
        self.locusdict = {}
        self.version = Version
        self.startenddict = {}
        self.transcriptnumlist = []
        self.seqdict = {}
        self.blastdict = {}
        self.compdict = {}

    def GFFparse(self):
        """This function parses through a GFF3 file to collect data from the genes, mRNA transcripts, and exons. If a gene,
        it collects gene name, gene ID, start and end of the gene, type of gene, and chromosome. If it were an mRNA transcript,
        the parser collected the gene ID, start and end of the transcript, transcript ID. Additionally, the number of iterations
        of the gene ID when parsing through mRNA transcripts was used to determine the number of transcripts for each gene.
        If it were an exon, the parser would collect the transcript ID and the number of iterations of the transcript ID was
        used to determine the number of exons for each transcript."""

        for i in self.GFF_file:
            if not i.startswith('#'):
                # get each line of information from the tab delimited file into a list for analysis
                line_list = i.strip().split('\t')
                if len(line_list) > 6:
                    # Get the information for each gene
                    if line_list[2] == 'CDS':
                        # Collect the start & end of the gene
                        Start = line_list[3]
                        End = line_list[4]
                        # get the attributes which have remaining info
                        attributes_list = line_list[8].split(';')
                        # get the gene name
                        Refseq = '0'
                        for i in attributes_list:
                            if i.startswith("inference"):
                                for j in i.split(':'):
                                    if j.startswith("WP"):
                                        Refseq = j
                                    # if j.startswith("WP") not in i.split(':'):
                                    # Refseq = "0"
                        Locus_list = attributes_list[1].split("-")
                        Locus = Locus_list[1]
                        # get the gene ID
                        ProteinID_list = attributes_list[0].split('=')
                        ProteinID = ProteinID_list[1].split('-')
                        ProteinID = ProteinID[1]
                        # get the gene product
                        for i in attributes_list:
                            if i.startswith("product="):
                                Product_list = i.split('=')
                                Product = Product_list[1]
                        self.locusdict[Locus] = (ProteinID, Locus, Start, End, Product, Refseq)
                        self.locusdict['0'] = 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'
    def FASTAparse(self):
        """This function parses through the FASTA to get the locus tag, protein name, start, end, gbkey, and
        sequence for each corresponding locus"""
        count = 0
        for i in self.FASTA_file:
                if  i.startswith('>'):
                    # get each line of information from the tab delimited file into a list for analysis
                    line_list = i.strip().split(' [')
                    gene=()
                    chromosome = ()
                    for j in line_list:
                        if j.startswith(">"):
                            blast_name = j.split('>')[1]
                        if j.startswith("locus_tag"):
                            locus = j.split('=')[1].strip(']')
                            self.seqdict[locus] = []
                        elif j.startswith('>'):
                            chromosome = j.split('|')[1].split('_')[0]
                        elif j.startswith("gene"):
                            gene = j.split('=')[1].strip(']')
                        elif j.startswith('protein='):
                            product = j.split('=')[1].strip(']')
                        elif j.startswith('location'):
                            start = j.split('=')[1].split('..')[0]
                            if start.startswith('complement('):
                                start = start.split('(')[1]
                            end = j.split('=')[1].split('..')[1].strip(']').strip(')')
                else:
                    addseq = i.strip()
                    self.seqdict[locus] = str(str(self.seqdict[locus]).strip('[').strip(']') + str(addseq))
                self.locusdict[locus] = locus, gene, product, start, end, self.seqdict[locus], chromosome
                self.locusdict['0'] = 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'
                self.blastdict[blast_name]=locus

    def analysisdf(self, dataframe):
        """This function takes a dataframe with at least one column labeled "Locus" and adds columns including the
        gene (if applicable), name, alternative locus, Start, End, Product, and chromosome from self.locusdict & a dictionary
        for A1S and ACX60 if chosen."""
        # add a column for each locus
        Locus_list = dataframe['Locus'].tolist()
        Gene_list = []
        Product_list = []
        Startend_list = []
        #chromosome_list = []
        for i in Locus_list:
            if i.startswith('gene-'):
                i = '_'.join(i.split('-')[1].split('_')[0:2])
            if i not in self.locusdict: i='0'
            Gene_list = Gene_list + [self.locusdict[i][1]]
            Product_list = Product_list + [self.locusdict[i][2]]
            Startend_list = Startend_list + [str(self.locusdict[i][3])+'...'+str(self.locusdict[i][4])]
            #chromosome_list = chromosome_list + [self.locusdict[i][6]]
        dataframe['Gene'] = Gene_list
        dataframe['Product'] = Product_list
        dataframe['Start, End'] = Startend_list
        #dataframe['Chrom'] = chromosome_list
        return dataframe



def genecomp(list, dict):
    """This function parses through a list and 'translates' it based on a dictionary, placing the translated form
    into a new list"""
    new_list = []
    for i in list:
        print('i', i)
        new_list = new_list + [dict[i]]
    return new_list

def genecompdf(df, dict):
    """This function takes a text delimited file with the A1S or ACX60 loci and adds a column with the opposite loci. It
     a new dataframe with the new column."""
    OrigLocus_list = df['Locus'].tolist()
    NewLocus_list = []
    for i in OrigLocus_list:
        if i in dict:
            NewLocus_list = NewLocus_list + [dict[i]]
        elif i.startswith('gene-A'):
            newi = '_'.join(i.split('-')[1].split('_')[0:2])
            if newi not in dict:
                NewLocus_list = NewLocus_list + ['NA']
            else:
                NewLocus_list = NewLocus_list + [dict[newi]]
        else: NewLocus_list = NewLocus_list + ['NA']
    df['New Locus'] = NewLocus_list
    return df

def compHB(list1, list2, only1, only2):
    """This function parses through two lists and determines which are shared and which are only found in list1
    or list2 and then returns two new lists corresponding to whether or not each observation in each list is
    "Shared" in both or only found in list1 (only1 as label) or in list2 (only2 as label)."""
    H_newlist = []
    B_newlist = []
    HBrel_dict = {}
    for i in list1:
        if i in list2:
            HBrel_dict[i] = 'Shared'
        else:
            HBrel_dict[i] = only1
    for j in list2:
        if j in list1:
            HBrel_dict[j] = 'Shared'
        else:
            HBrel_dict[j] = only2
    for m in list1:
        H_newlist = H_newlist + [HBrel_dict[m]]
    for n in list2:
        B_newlist = B_newlist + [HBrel_dict[n]]
    return H_newlist, B_newlist

def Blastcomp1(Blast_file):
    '''This function takes 1 blast file (txt) and returns a dictionary comparing the Query & Subject'''
    Blastd = pd.read_csv(Blast_file, sep="\t",
                         names=['Query', 'Subject', '% alignment', 'length', '# of mismatches',
                                '# of gaps', 'qstart', 'qend', 'sstart', 'send', 'E-value', 'bitscore'])
    Blastdf = pd.DataFrame(Blastd)
    Query_list = Blastd['Query'].tolist()
    Subject_list = Blastd['Subject'].tolist()
    blastdict = {}

    for (i, j) in zip(Query_list, Subject_list):
        if i.startswith('lcl|NZ'): i = '_'.join(i.split('_')[3:5])
        else: i = i.split('_')[2]
        if j.startswith('lcl|NZ'): j = '_'.join(j.split('_')[3:5])
        else: j = j.split('_')[2]
        blastdict[i] = j
    return blastdict

def get_key(val, my_dict):
    for key, value in my_dict.items():
        if val == value:
            return key



if __name__ == "__main__":

    mffrs_a1s_dict = Blastcomp1('Genomes/blasta1smff.txt')
    mffrs_mff_dict = Blastcomp1('Genomes/BlastRefSeqMff.txt')
    A1S = Genome(GFF3file="Genomes/A1S.gff3", Fastafile="Genomes/A1S.txt", Version='a1s')
    A1S.GFFparse()


    Mff = Genome(GFF3file="Genomes/mfffull.gff3", Fastafile="Genomes/mfffullfasta.txt", Version='mff')
    Mff.GFFparse()


    newtestdict = {}
    newtestdict['0'] = 'NA'
    for key, value in mffrs_mff_dict.items():
        print('key', key, 'value', value)
         for i in A1S.locusdict.values():
             proteinid = i[0] #(ABO or AKQ)
             locus = i[1]
             if key in newtestdict.keys():
                 newtestdict[key] = newtestdict[key] + [locus, value]
             else:
                 newtestdict[key] = [locus, value]
        for j in Mff.locusdict.values():
            proteinid1 = j[0]
            locus1 = j[1]
            if proteinid1 == value:
                newtestdict[locus1] = key
            if key in newtestdict.keys():
                newtestdict[key] = newtestdict[key] + [locus1, value]


    for key, value in mffrs_a1s_dict.items():
        for l in A1S.locusdict.values():
            #print(i)
            proteinid = l[0] #(ABO or AKQ)
            #print(proteinid)
            locus2 = l[1]
            if proteinid == value:
                if key in newtestdict.keys():
                    newtestdict[key] = newtestdict[key] + [locus2, value]
                else: newtestdict[key] = [locus2, value]
        for m in Mff.locusdict.values():
            proteinid1 = m[0]
            #print(proteinid1)
            locus3 = m[1]
            if proteinid1 == value:
                if key in newtestdict.keys():
                    newtestdict[key] = newtestdict[key] + [locus3, value]
                else:
                    newtestdict[key] = [locus3, value]



    """Get the dictionary to compare A1S & ACX60 genomes"""
    Blastd = pd.read_csv("Genomes/blasta1smff.txt", sep="\t", names=['Query','Subject', '% alignment', 'length', '# of mismatches',
                                                                '# of gaps', 'qstart', 'qend', 'sstart', 'send', 'E-value', 'bitscore'])


    comparedict = {}
    for (i, j) in zip(NewQuery_list, NewSubject_list):
        comparedict[j]= i


    """Create a massive csv file with the appropriate info for ∆BlsA and WT csv files from R"""

    # Create dataframes from the files for BlsA & WT
    BlsAd = pd.read_csv("RNASeqPhrB_files/Part2/FinalMatrixDEBlsA2.csv", sep=",", header=0,
                    names=['Number', 'Locus', 'PPEE', 'PPDE', 'PostFC', 'RealFC',
                          'Direction', 'D1', 'D2', 'D3', 'L1', 'L2', 'L3'])
    BlsAdf = pd.DataFrame(BlsAd)
    # iterate through the Locus file and make a new dataframe with an additional column
    BlsAnewdf = BlsAdf
    # BlsAnewdf = genecompdf(BlsAdf, comparedict)
    B_list = BlsAnewdf['Locus'].tolist()
    Andres_list = []
    for i in B_list:
        i = '_'.join(i.split('-')[1].split('_')[:2])
        if i not in newtestdict: i= '0'
        Andres_list = Andres_list + [newtestdict[i]]
    BlsAnewdf['Andres'] = Andres_list


    WTd = pd.read_csv("RNASeqPhrB_files/Part2/FinalMatrixDEWT2.csv", sep=",",
                     names=['Number', 'Locus', 'PPEE', 'PPDE', 'PostFC', 'RealFC', 'Direction',
                            'D1', 'D2', 'D3', 'L1', 'L2', 'L3'])
    WTdf = pd.DataFrame(WTd)
    WTnewdf = WTdf
    WTnewdf = genecompdf(WTdf, comparedict)
    WT_list = WTnewdf['Locus'].tolist()


    # create lists for the genes in BlsA and WT
    WT_newlist, BlsA_newlist = compHB(WT_list, B_list, 'WT only', 'BlsA only')
    # Add a column to state whether or not the gene is only DE in BlsA or both
    WTnewdf['Relation'] = WT_newlist
    BlsAnewdf['Relation'] = BlsA_newlist
    # Add the RefSeq one onto Hope's data
    RefseqW_list = [0]
    RefseqB_list = []
    for i in WT_list[1:]:
        i = '_'.join(i.split('-')[1].split('_')[:2])
        if i in Mff.locusdict:
            RefseqW_list = RefseqW_list + [Mff.locusdict[i][5]]
        else: RefseqW_list = RefseqW_list + ['0']
    for i in B_list:
        i = '_'.join(i.split('-')[1].split('_')[:2])
        if i in Mff.locusdict:
            RefseqB_list = RefseqB_list + [Mff.locusdict[i][5]]
        else: RefseqB_list = RefseqB_list + ['0']
    # in Blsa --> 'BlsA only'

    # Add columns to include other things
    BlsAnewdf = Mff.analysisdf(dataframe=BlsAnewdf)
    BlsAnewdf['RefSeq'] = RefseqB_list
    WTnewdf = Mff.analysisdf(dataframe=WTnewdf)
    WTnewdf['RefSeq'] = RefseqW_list


    # Add a column to state whether or not the gene is only DE in BlsA or both
    WTnewdf['H/B Relation'] = WTH_newlist
    BlsAnewdf['H/B Relation'] = BH_newlist
    # save file
    WTnewdf.to_csv('RNASeqPhrB_files/Part2/FinalWT.csv')
    BlsAnewdf.to_csv('RNASeqPhrB_files/Part2/FinalBlsA.csv')
  





