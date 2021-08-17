import pandas as pd


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
    vep_dataframe = vep_dataframe[vep_dataframe['%SYMBOL_SOURCE'] == 'HGNC']
    vep_dataframe.drop_duplicates(inplace=True)
    vep_dataframe['chr_pos_ref_alt'] = vep_dataframe[[
        '%CHROM', '%POS', '%REF', '%ALT'
    ]].astype(str).agg('_'.join, axis=1)
    variant_dataframe['chr_pos_ref_alt'] = variant_dataframe[[
        '#CHROM', 'POS', 'REF', 'ALT'
    ]].astype(str).agg('_'.join, axis=1)
    merge = vep_dataframe.merge(variant_dataframe, on='chr_pos_ref_alt')
    merge.drop(
        columns=['chr_pos_ref_alt', '#CHROM', 'POS', 'REF', 'ALT'],
        inplace=True
    )
    merge['binarized_label'] = 0
    merge.loc[merge[merge['clinsig'] == 'LP'].index, 'binarized_label'] = 1
    return merge
