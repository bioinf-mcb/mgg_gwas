import pandas as pd
import numpy as np
import re

def report_pfam(df, verbose=True):
    
    # get PC name
    pcid = df['query'].unique()[0]
    
    # sort hits
    informative_df = df.query('~name.str.contains("DUF")').sort_values('bits', ascending=False)
    noninformative_df = df.query('name.str.contains("DUF")').sort_values('bits', ascending=False)
    df = pd.concat([informative_df, noninformative_df])
    
    # best hit & function2report
    if len(df) != 0: 
        top_hit_df = df.iloc[0].copy()
        name = top_hit_df['name']
        pfamID, func_short, func_detailed = name.split(';')[0].strip(), name.split(';')[1].strip(), name.split(';')[2].strip()
        report_function = f'{func_detailed} [{func_short}] [{pfamID}]'
        report_params = f'bits: {int(top_hit_df["bits"]):>5}  eval: {top_hit_df["evalue"]:>6.1E}'
    else: # no hit 
        top_hit_df = get_no_hit_row(1).iloc[0]
        report_function, report_params = '-', '-'
    
    
    # report columns
    top_hit_df['query'] = pcid
    top_hit_df['report_label'] = 'PFAM'
    top_hit_df['report_function'] = report_function
    top_hit_df['report_params'] = report_params
    top_hit_df['report_confidence'] = get_confidence_column(top_hit_df)
            
    if verbose: display(top_hit_df.to_frame().T)
    return top_hit_df.to_frame().T
    

def report_unique_ecod_topologies(ecod_df, pcs80):
    """ for each PC report 5 unique topologies for different X levels with max bitscore """


    # ECOD regex
    def parse_ecod_info(row, pattern_ecodID = r"\|\s*((?:\d+\.\d+\.\d+\.\d+)|(?:\d+\.\d+\.\d+))\s*\|"):
        
        # default
        ecodID, A_LEVEL, X_LEVEL, H_LEVEL, T_LEVEL, F_LEVEL, PROTEIN_LEVEL = ['EXTRACTION_ERROR'] * 7


        try:
            # initial split
            ecod_id_raw, rest = row['name'].split(' | A: ')

            # ECOD ID
            match = re.search(pattern_ecodID, ecod_id_raw)
            ecodID = match.group(1)

            # ECOD levels
            A_LEVEL, rest = rest.split(', X: ')
            X_LEVEL, rest = rest.split(', H: ')
            H_LEVEL, rest = rest.split(', T: ')
            T_LEVEL, rest = rest.split(', F: ')
            F_LEVEL, PROTEIN_LEVEL = rest.split(' | ')
            PROTEIN_LEVEL = PROTEIN_LEVEL.strip('Protein: ')

        except: pass
        
        return pd.Series([ecodID, f'A: {A_LEVEL}', f'X: {X_LEVEL}', f"H: {H_LEVEL}", f"T: {T_LEVEL}", f"F: {F_LEVEL}", f"PROTEIN: {PROTEIN_LEVEL}"])

    # apply
    parsed_ecod = ecod_df.apply(parse_ecod_info, axis=1)
    parsed_ecod.columns = ['ecodID', 'A_LEVEL', 'X_LEVEL', 'H_LEVEL', 'T_LEVEL', 'F_LEVEL', 'PROTEIN_LEVEL']

    # combine
    ecod_df = pd.concat([ecod_df, parsed_ecod], axis=1)

    # clean ecod
    ecod_cols = ['PC80', 'X_LEVEL', 'T_LEVEL', 'F_LEVEL', 'qstart', 'qend', 'prob', 'bits', 'ecodID']
    ecod_df = ecod_df.rename({'query': 'PC80'}, axis=1)
    ecod_df['PC80-sort'] = pd.to_numeric(ecod_df['PC80'].str.strip('PC'), downcast='integer')
    ecod_df = ecod_df.sort_values(['PC80-sort', 'bits'], ascending=[True, False])
    ecod_df = ecod_df[ecod_cols]


    ### Report unique topologies for PC80        

    def get_unique_t_levels(group):
        """ report 5 unique topologies per PC80 """

        # sort & select
        group = group.sort_values('prob', ascending=False)
        idx = group.groupby('X_LEVEL')['bits'].idxmax()
        unique_t_levels = group.loc[idx, 'T_LEVEL'].to_list()

        # holder
        unique_t_levels = unique_t_levels + ['-'] * 5
        
        return unique_t_levels[:5]

    unique_topologies = []
    for pcid, group in ecod_df.groupby('PC80'):
        topology1, topology2, topology3, topology4, topology5 = get_unique_t_levels(group)
        row = {'PC80': pcid, 'ECOD1': topology1,'ECOD2': topology2, 'ECOD3': topology3, 'ECOD4': topology4, 'ECOD5': topology5}
        unique_topologies.append(row)

    # add missing pcs
    pcs80 = [pc for pc in pcs80 if pc not in ecod_df['PC80'].unique()]
    for pcid in pcs80:
        topology1, topology2, topology3, topology4, topology5 = ['-'] * 5
        row = {'PC80': pcid, 'ECOD1': topology1,'ECOD2': topology2, 'ECOD3': topology3, 'ECOD4': topology4, 'ECOD5': topology5}
        unique_topologies.append(row)

    # create
    unique_topologies_df = pd.DataFrame(unique_topologies)

    # sort
    unique_topologies_df['PC80-sort'] = pd.to_numeric(unique_topologies_df['PC80'].str.strip('PC'), downcast='integer')
    unique_topologies_df = unique_topologies_df.sort_values('PC80-sort', ascending=True)
    unique_topologies_df = unique_topologies_df.drop('PC80-sort', axis=1)

    return unique_topologies_df

def get_unique_functions_frame(df, function_column='annot'):
    """ for each function in data frame:
    - get best hit [highest bitscore]
    - for get highest qcov for a given function from all of the hits
    - return frame of best hits for each function with highest qcov """
    
    # select best hits for each unique function (highest bitscore)
    best_hits_df = df.loc[df.groupby(function_column)['bits'].idxmax()] \
                                                             .sort_values('bits', ascending=False) \
                                                             .copy()
    
    # get highest qcov for each unique function
    best_qcov_df = df.loc[df.groupby(function_column)['qcov'] \
                            .idxmax()][[function_column,'qcov']] \
                            .copy()
    
    # best hits for each unique function & highest qcov for given function
    final_df = best_hits_df.merge(best_qcov_df, on=function_column, how='left', suffixes=('_oryginal', '_best')) \
                           .drop('qcov_oryginal', axis=1) \
                           .rename(columns={'qcov_best': 'qcov'}) \
                           .sort_values('bits', ascending=False)
    
    return final_df

def report_phrogs(df, max_evalue=10**-3, nfunc2report=2, verbose=False):

    """
    report a specific number (nfunc2report) of unique functions that got significat (max_evalue) hit(s).
    However, firsty report function significant (not in non-informative list of functions), secondly non-informative functions, lastly unknown functions.
    
    Algorithm:
    1. Filter eval 10**-3
    2. Remove unknown function [REPORT ONLY WHEN NO OTHER FUNCTION, PRIORITY-0]
    3. Remove non-informative functions (lytic tail protein, tail protein, structural protein, virion structural protein, minor tail protein ...) [REPORT ONLY WHEN NO OTHER FUNCTION; PRIORITY-1]
    4. Group by unique functions. For each function report independently max bitscore and max qcov (hits to this function).
    5. Take two functions with highest bitscores.
    6. Report {confidence} {function} in one genbank field, and in seperate field bitscore and qcov.
    7. Report two best PHROG hits seperataly (in total four PHROGS field: 2x function with confidence and 2x params: bitscore and qcov) [PRIORITY-2]
    
    """
    
    # get PC name
    pcid = df['query'].unique()[0]
    
    ### get filters
    noninformative_functions = ['lytic tail protein', 'tail protein', 'structural protein', 'virion structural protein', 'minor tail protein']

    filt_evalue = 'evalue <= @max_evalue'
    get_unknown = '(annot == "unknown function")'
    get_noninformative_functions = '(annot.isin(@noninformative_functions))'
    remove_unknown = '~' + get_unknown
    remove_noninformative = '~' + get_noninformative_functions
    
    informative_query = ' and '.join([remove_unknown, remove_noninformative])
    noninformative_query = get_noninformative_functions

    ### significat only
    df = df.query(filt_evalue)
    
    ### get informative hits
    informative_df = df.query(informative_query)
    informative_df = get_unique_functions_frame(informative_df, function_column='annot') # best hit [max bit score] & highest qcov

    ### noninformative hits
    noninformative_df = df.query(get_noninformative_functions)
    noninformative_df = get_unique_functions_frame(noninformative_df, function_column='annot') # best hit [max bit score] & highest qcov
        
    ### uknown hits
    unknown_df = df.query(get_unknown)
    unknown_df = get_unique_functions_frame(unknown_df, function_column='annot') # best hit [max bit score] & highest qcov

    ### report function
    top_hits_df = pd.concat([informative_df, noninformative_df, unknown_df]).iloc[:nfunc2report].copy()
    
    # prepare columns2report
    if len(top_hits_df) == 0:
        top_hits_df = get_no_hit_row(2)    
        top_hits_df['report_params'] = ['-', '-']
    elif len(top_hits_df) == 1: 
        holder_row = get_no_hit_row(1)
        holder_row['report_params'] = '-'

        bits, qcov = top_hits_df["bits"].iloc[0], top_hits_df["qcov"].iloc[0]
        top_hits_df['report_params'] = f'bits: {int(bits):>5}  qcov: {qcov:.2f}'
        top_hits_df = pd.concat([top_hits_df, holder_row])
    else: 
        top_hits_df['report_params'] = top_hits_df.apply(lambda row: f'bits: {int(row["bits"]):>5}  qcov: {row["qcov"]:>6.2f}', axis=1)

    # report columns
    labels = [f'PHROGS{i}' for i in range(1,len(top_hits_df)+1)]
    
    top_hits_df['query'] = [pcid] * len(top_hits_df)
    top_hits_df['report_label'] = labels
    top_hits_df['report_function'] = top_hits_df['annot']
    top_hits_df['report_confidence'] = get_confidence_column(top_hits_df)
            
    if verbose: display(top_hits_df)
    return top_hits_df


def get_confidence_column(df, eval_intervals=(10**-10, 10**-5, 10**-3, 1), col='evalue', verbose=False):
    # Build conditions: first for evalue == 0, then for each interval, and one for values >= the last interval.
    conditions = []
    
    # Condition for evalue equal to 0: assign '***'
    conditions.append(df[col] == 0)
    
    # For the given intervals
    for i, threshold in enumerate(eval_intervals):
        if i == 0:
            conditions.append(df[col] <= threshold)
        else:
            lower = eval_intervals[i-1]
            conditions.append((df[col] > lower) & (df[col] <= threshold))
    
    # Condition for evalue above the last threshold.
    conditions.append(df[col] >= eval_intervals[-1])
    
    # Generate star symbols for the intervals.
    # For example, with 4 intervals, star_choices becomes: ['***', '**', '*', ''].
    star_choices = ['*' * i for i in range(len(eval_intervals)).__reversed__()]
    
    # Now, create the final list of choices:
    # - For evalue==0: '***'
    # - For the best interval (lowest evalue): use the first star value from star_choices (if any)
    # - Then the rest of the intervals, and finally two conditions for high evalues.
    choices = ['***'] + star_choices[:-1] + ['!', '!']
    
    if verbose:
        print("Conditions:", conditions)
        print("Choices:", choices)
    
    return np.select(conditions, choices, default='?')


def get_no_hit_row(n, columns_mapper = {'query': 'string', 'target': 'string', 'prob': 'float', \
                                        'pvalue': 'float', 'ident': 'float', 'qcov': 'float', \
                                        'tcov': 'float', 'bits': 'float', 'qstart': 'int', \
                                        'qend': 'int', 'qlength': 'int', 'tstart': 'int', \
                                        'tend': 'int', 'tlength': 'int', 'evalue': 'float', \
                                        'db': 'string', 'name': 'string', 'color': 'string', \
                                        'annot': 'string', 'category': 'string', 'phrog/alan_profile': 'string',
                                        'report_label': 'string', 'report_function': 'string', 'report_params': 'string'}):
    
    """ Give dict of column names and variable types to create 'no hit' row as data frame object """

    values, indicies = [], columns_mapper.keys()
    for key, variable_type in columns_mapper.items():
        if variable_type == 'string': values.append('-')
        else: values.append(0)

    no_hit_row = pd.Series(values, index=indicies).to_frame().T
    no_hit_row = pd.concat([no_hit_row]*n)
    return no_hit_row


