import pandas as pd
import gzip
import os


def order_vcf(vcf):
    """
    Function to order a TSV read in VCF file by chromosome 1 to 22, X, Y and
    MT.

    Parameters
    ----------
    vcf : pandas.DataFrame
        The as TSV read in pandas dataframe of the VCF

    Returns
    -------
    pandas.DataFrame
        Ordered VCF in TSV format, ordered by 1 through 22, X, Y and MT.
    """
    vcf['order'] = vcf['#CHROM']
    vcf.loc[vcf[vcf['#CHROM'] == 'X'].index, 'order'] = 23
    vcf.loc[vcf[vcf['#CHROM'] == 'Y'].index, 'order'] = 24
    vcf.loc[vcf[vcf['#CHROM'] == 'MT'].index, 'order'] = 25
    vcf['order'] = vcf['order'].astype(int)
    vcf.sort_values(by=['order', 'POS'], inplace=True)
    vcf.drop(columns='order', inplace=True)
    vcf.index = range(0, vcf.shape[0])
    return vcf


def merge_vep_variant_files(vep_dataframe: pd.DataFrame,
                            variant_dataframe: pd.DataFrame):
    """
    Function to merge the train / validation dataframe output from VEP with
    the original variant csv.

    Parameters
    ----------
    vep_dataframe: pandas.DataFrame
        Loaded in pandas dataframe of the raw VEP output VCF.

    variant_dataframe: pandas.DataFrame
        Loaded in pandas dataframe of the variant csv file.

    Returns
    -------
    merge: pandas.DataFrame
        Merged pandas dataframe of the VEP output with the added "clinsig",
        "review", "stars" and "binarized_label" columns to the vep_dataframe.
    """
    print('Starting merge.')
    vep = vep_dataframe.copy(deep=True)
    vep.drop(
        index=vep[vep['%SYMBOL_SOURCE'] != 'HGNC'].index,
        inplace=True
    )
    print('Dropping duplicates.')
    vep.drop_duplicates(inplace=True)
    print('Creating merge column.')
    vep['chr_pos_ref_alt'] = vep[[
        '%CHROM', '%POS', '%REF', '%ALT'
    ]].astype(str).agg('_'.join, axis=1)
    variant_dataframe['chr_pos_ref_alt'] = variant_dataframe[[
        '#CHROM', 'POS', 'REF', 'ALT'
    ]].astype(str).agg('_'.join, axis=1)
    print('Merging.')
    merge = vep.merge(variant_dataframe, on='chr_pos_ref_alt')
    merge.drop(
        columns=['chr_pos_ref_alt', '#CHROM', 'POS', 'REF', 'ALT'],
        inplace=True
    )
    print('Merge done. Adding binarized label.')
    merge['binarized_label'] = 0
    merge.loc[merge[merge['clinsig'] == 'LP'].index, 'binarized_label'] = 1
    print('Done.')
    return merge


class SimpleVCFParser:
    def __init__(self, vcf_file: str):
        """
        SimpleVCFParser, for when you want to read in a VCF file quickly
        through pandas.

        Parameters
        ----------
        vcf_file: str, path-like
            Location to the gzipped VCF file.
        """
        self.infile = vcf_file
        self.header = []
        self.data = pd.DataFrame()
        if self.infile:
            self._check_tilde()
            self._read_header()
            self._read_data()

    def _check_tilde(self):
        if self.infile.startswith('~'):
            self.infile = os.path.expanduser(self.infile)

    def _read_header(self):
        conn = gzip.open(self.infile, 'rt')
        for line in conn:
            if line.startswith('##'):
                self.header.append(line)
            else:
                break
        conn.close()

    def _read_data(self):
        self.data = pd.read_csv(
            self.infile, sep='\t', skiprows=len(self.header), low_memory=False
        )

    @property
    def data(self):
        """
        Property to get the variant data within the VCF.

        Returns
        -------
        data: pandas.DataFrame
            Pandas dataframe of the variant data.
        """
        return self._data

    @data.setter
    def data(self, value):
        """
        Setter data, to overwrite the data currently stored in the class.

        Parameters
        ----------
        value: pandas.DataFrame
            Pandas dataframe of the variant data. Doesn't have to be VCF format.

        Returns
        -------
        data: pandas.DataFrame
            Pandas dataframe of the variant data.

        Raises
        ------
        ValueError:
            If value is not a pandas DataFrame.
        """
        if not isinstance(value, pd.DataFrame):
            raise ValueError('data can not be anything other than a pandas '
                             'dataframe!')
        self._data = value

    def export(self, file_loc):
        """
        Exporter function of SimpleVCFParser to write the header to a file
        followed by the property data. Can be either TSV format or VCF.

        Parameters
        ----------
        file_loc: str, path-like
            Location of where the output file should be. Will always be gzipped.

        Raises
        ------
        ValueError
            When data and/or header is not set.
        """
        if self._data.shape[0] == 0 or len(self.header) == 0:
            raise ValueError('data and header have to be set in order '
                             'to export anything!')
        if file_loc.startswith('~'):
            file_loc = os.path.expanduser(file_loc)
        if not os.path.exists(os.path.dirname(file_loc)):
            print('Unable to locate path, creating.')
            os.makedirs(os.path.dirname(file_loc))
            print('Output directory made.')
        if not file_loc.endswith('.gz'):
            file_loc += '.gz'
        with gzip.open(file_loc, 'wt') as output_file:
            for line in self.header:
                output_file.write(line)
        self.data.to_csv(file_loc, sep='\t', mode='a', index=False)
        print('Successfully exported to: {}'.format(file_loc))


pseudo_header = ['##fileformat=VCFv4.2',
                 '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">',
                 '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                 '##contig=<ID=1,length=249250621,assembly=b37>',
                 '##contig=<ID=2,assembly=b37,length=243199373>',
                 '##contig=<ID=3,assembly=b37,length=198022430>',
                 '##contig=<ID=4,length=191154276,assembly=b37>',
                 '##contig=<ID=5,length=180915260,assembly=b37>',
                 '##contig=<ID=6,length=171115067,assembly=b37>',
                 '##contig=<ID=7,length=159138663,assembly=b37>',
                 '##contig=<ID=8,length=146364022,assembly=b37>',
                 '##contig=<ID=9,length=141213431,assembly=b37>',
                 '##contig=<ID=10,length=135534747,assembly=b37>',
                 '##contig=<ID=11,length=135006516,assembly=b37>',
                 '##contig=<ID=12,length=133851895,assembly=b37>',
                 '##contig=<ID=13,length=115169878,assembly=b37>',
                 '##contig=<ID=14,length=107349540,assembly=b37>',
                 '##contig=<ID=15,length=102531392,assembly=b37>',
                 '##contig=<ID=16,length=90354753,assembly=b37>',
                 '##contig=<ID=17,length=81195210,assembly=b37>',
                 '##contig=<ID=18,length=78077248,assembly=b37>',
                 '##contig=<ID=19,length=59128983,assembly=b37>',
                 '##contig=<ID=20,length=63025520,assembly=b37>',
                 '##contig=<ID=21,length=48129895,assembly=b37>',
                 '##contig=<ID=22,length=51304566,assembly=b37>',
                 '##contig=<ID=X,assembly=b37,length=155270560>',
                 '##contig=<ID=Y,length=59373566,assembly=b37>',
                 '##contig=<ID=MT,length=16569,assembly=b37>',
                 '##fileDate=20200320',
                 '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|UID">']
