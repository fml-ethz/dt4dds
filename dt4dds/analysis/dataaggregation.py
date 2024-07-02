import numpy as np
import pandas as pd
import pathlib
import ruamel.yaml
import Bio.SeqIO
import Bio.SeqUtils

rng = np.random.default_rng()


def read_reference(input_file):
    d = dict()
    with open(input_file) as f:
        for record in Bio.SeqIO.parse(f, 'fasta'):
            d[str(record.id)] = record
    return d



class ErrorFile():

    def __init__(self, file, read=1, skip=[]):

        self.file = pathlib.Path(file)
        self.read = read
        yaml = ruamel.yaml.YAML(typ='safe', pure=True)
        data = yaml.load(open(file))
        self.raw_data = {}
        self.data = {}
        self.exp = self.file.parent.parent.name

        for var, val in data.items():

            if type(val) == dict:
                dict_contents = list(val.items())

                if len(dict_contents) > 0 and type(dict_contents[0][1]) == dict:
                    formatted_val = {t: pd.DataFrame.from_dict(v, orient='index', columns=['value']) for t,v in val.items()}
                else:
                    formatted_val = pd.DataFrame.from_dict(val, orient='index', columns=['value'])

            else:
                formatted_val = val


            setattr(self, var, formatted_val)
            self.raw_data[var] = formatted_val

        
        # PROPERTIES
        self.data['n_reads'] = pd.DataFrame.from_dict({'n_reads': [self.raw_data['n_reads']]}, orient='columns')
        
        # OVERALL ERROR RATES
        df = pd.DataFrame({
            'type': ('substitutions', 'deletions', 'insertions', 'delevents', 'subevents', 'insevents'),
            'rate': [getattr(self, f"r_{errortype}") for errortype in ('substitutions', 'deletions', 'insertions', 'delevents', 'subevents', 'insevents')]
        })
        df['read'] = str(self.read)
        df['exp'] = self.exp
        self.data[f"overall_error_rates"] = df

        # BASES BY POSITION
        for basetype in ('n_bases_by_position', 'n_bases_by_refposition'):
            df = self.raw_data[basetype].copy().reset_index()
            df = df.rename(columns={'index': 'position', 'value': 'bases'})
            df['read'] = str(self.read)
            df['exp'] = self.exp
            if skip:
                df.drop(df[(df['position'] >= skip[0]) & (df['position'] <= skip[1])].index, inplace=True)
            self.data[basetype] = df

        # OVERALL ERROR BIASES
        for errortype in ('substitutions', 'deletions', 'insertions'):
            df = self.raw_data[f"p_{errortype}_by_type"].copy().reset_index()
            df = df.rename(columns={'index': 'type', 'value': 'ratio'})
            df['rate'] = df['ratio']*self.raw_data[f"r_{errortype}"]
            df['read'] = str(self.read)
            df['exp'] = self.exp
            self.data[f"{errortype}_by_type"] = df

        # ERROR RATES BY POSITION
        for errortype in ('substitutions', 'deletions', 'insertions'):
            df = self.raw_data[f"p_{errortype}_by_position"].copy().reset_index()
            df = df.rename(columns={'index': 'position', 'value': 'rate'})
            df = df.drop(df[df['position'] < 0].index)
            df['read'] = str(self.read)
            df['exp'] = self.exp
            if skip:
                df.drop(df[(df['position'] >= skip[0]) & (df['position'] <= skip[1])].index, inplace=True)
            self.data[f"{errortype}_by_position"] = df

        # ERROR RATES BY REFPOSITION
        for errortype in ('substitutions', 'deletions', 'insertions'):
            df = self.raw_data[f"p_{errortype}_by_refposition"].copy().reset_index()
            df = df.rename(columns={'index': 'position', 'value': 'rate'})
            df = df.drop(df[df['position'] < 0].index)
            df['read'] = str(self.read)
            df['exp'] = self.exp
            if skip:
                df.drop(df[(df['position'] >= skip[0]) & (df['position'] <= skip[1])].index, inplace=True)
            self.data[f"{errortype}_by_refposition"] = df

        # COMBINED ERROR RATES BY POSITION
        dfs = []
        for errortype in ('substitutions', 'deletions', 'insertions'):
            df = self.data[f"{errortype}_by_position"].copy()
            df['type'] = errortype
            if not df.empty: dfs.append(df)
        self.data[f"overall_error_rates_by_position"] = pd.concat(dfs).reset_index(drop=True)

        # COMBINED ERROR RATES BY REFPOSITION
        dfs = []
        for errortype in ('substitutions', 'deletions', 'insertions'):
            df = self.data[f"{errortype}_by_refposition"].copy()
            df['type'] = errortype
            if not df.empty: dfs.append(df)
        self.data[f"overall_error_rates_by_refposition"] = pd.concat(dfs).reset_index(drop=True)

        # ERROR BIAS BY POSITION
        for errortype in ('substitutions', 'deletions', 'insertions'):
            dfs = []
            data = self.raw_data[f"p_{errortype}_by_position_by_type"]
            for t, df in data.items():
                df = df.copy().reset_index()
                df = df.rename(columns={'index': 'position', 'value': 'rate'})
                df = df.drop(df[df.position < 0].index)
                df['type'] = t
                df['read'] = str(self.read)
                df['exp'] = self.exp
                if skip:
                    df.drop(df[(df['position'] >= skip[0]) & (df['position'] <= skip[1])].index, inplace=True)
                if not df.empty: dfs.append(df)
            all_data = pd.concat(dfs).reset_index(drop=True)
            all_data['ratio'] = all_data['rate']/all_data.groupby('position')['rate'].transform("sum")
            self.data[f"{errortype}_by_position_by_type"] = all_data

        # ERROR BIAS BY REFPOSITION
        for errortype in ('substitutions', 'deletions', 'insertions'):
            dfs = []
            data = self.raw_data[f"p_{errortype}_by_refposition_by_type"]
            for t, df in data.items():
                df = df.copy().reset_index()
                df = df.rename(columns={'index': 'position', 'value': 'rate'})
                df = df.drop(df[df.position < 0].index)
                df['type'] = t
                df['read'] = str(self.read)
                df['exp'] = self.exp
                if skip:
                    df.drop(df[(df['position'] >= skip[0]) & (df['position'] <= skip[1])].index, inplace=True)
                if not df.empty: dfs.append(df)
            all_data = pd.concat(dfs).reset_index(drop=True)
            all_data['ratio'] = all_data['rate']/all_data.groupby('position')['rate'].transform("sum")
            self.data[f"{errortype}_by_refposition_by_type"] = all_data

        # ERROR FREQUENCY BY LENGTH
        dfs = []
        for errortype in ('substitutions', 'deletions', 'insertions'):
            df = self.raw_data[f"p_{errortype}_by_length"].reset_index()
            df = df.rename(columns={'index': 'length'})
            df['type'] = errortype
            df['read'] = str(self.read)
            df['exp'] = self.exp
            if not df.empty: dfs.append(df)
        self.data[f"error_frequency_by_length"] = pd.concat(dfs).reset_index(drop=True)

        # ERROR FREQUENCY PER READ
        dfs = []
        for errortype in ('substitutions', 'deletions', 'insertions', 'delevents', 'subevents', 'insevents'):
            df = self.raw_data[f"p_{errortype}_by_read"].reset_index()
            df = df.rename(columns={'index': 'frequency'})
            df['type'] = errortype
            df['read'] = str(self.read)
            df['exp'] = self.exp
            if not df.empty: dfs.append(df)
        self.data[f"error_frequency_by_read"] = pd.concat(dfs).reset_index(drop=True)






class ErrorAnalysis():
    
    def __init__(self, folder, skip=[(), ()], local=False, paired=False, category="mapped_high"):

        # build the file path to the analysis files
        r1_path = f"fw.{'local' if local else 'global'}.{category}.stats"
        r2_path = f"rv.{'local' if local else 'global'}.{category}.stats"
        self.folder = pathlib.Path(folder)

        # read analysis files
        reads = []
        self.r1 = ErrorFile(self.folder / r1_path, read=1, skip=skip[0])
        reads.append(self.r1)
        if paired:
            self.r2 = ErrorFile(self.folder / r2_path, read=2, skip=skip[1])
            reads.append(self.r2)

        # combine the data from both reads
        self.variables = list(self.r1.data.keys())
        self.data = {}
        for var in self.variables:
            dfs = [r.data[var] for r in reads]
            self.data[var] = pd.concat(dfs).reset_index(drop=True)




class GroupAnalysis():

    def __init__(self, experiments):

        self.experiments = experiments

        # combine the data from all analysis objects and add group identifier
        self.variables = experiments[0][-1].variables
        self.data = {}
        for var in self.variables:
            dfs = []
            for group, erroranalysis in self.experiments:
                df = erroranalysis.data[var].copy()
                df['group'] = group
                dfs.append(df)
            self.data[var] = pd.concat(dfs).reset_index(drop=True)





class SeriesAnalysis():


    def __init__(self, experiments):

        self.experiments = experiments

        # combine the data from all analysis objects and add group identifier
        self.variables = experiments[0][-1].variables
        self.data = {}
        for var in self.variables:
            dfs = []
            for group, series_var, erroranalysis in self.experiments:
                df = erroranalysis.data[var].copy()
                df['group'] = group
                df['series_var'] = series_var
                dfs.append(df)
            data = pd.concat(dfs).reset_index(drop=True)

            data['delta_series_var'] = data['series_var']
            for group in data['group'].unique():
                data.loc[(data['group'] == group), 'delta_series_var'] -= min(data.loc[(data['group'] == group), 'series_var'])
            
            self.data[var] = data



class SeqGroupAnalysis():

    def __init__(self, experiments):
        
        self.experiments = experiments

        dfs = []
        for group, erroranalysis in self.experiments:
            df = erroranalysis.overview_data.copy()
            df['group'] = group
            dfs.append(df)
        self.overview_data = pd.concat(dfs).reset_index(drop=True)

        self.position_data = {}
        for group, erroranalysis in self.experiments:
            if group not in self.position_data:
                self.position_data[group] = []
            self.position_data[group].append(erroranalysis)




class DistributionAnalysis():

    def __init__(self, samples):

        self.samples = samples
        datas = []

        for exp, file in samples.items():
            data = pd.read_csv(file, sep="\t", converters={'#name': str})[['#name', 'assignedReads']]
            refset = set(read_reference(pathlib.Path(file).parent / "../design_files.fasta").keys())
            missing = set(data['#name']) ^ refset
            d = pd.DataFrame.from_dict({
                '#name': list(missing), 
                'assignedReads': np.zeros(len(missing)),
                'x': np.zeros(len(missing)),
            })
            print(f"Missing sequences: {len(missing)}, {len(missing)/len(refset)*100:.2f}%")
            data['x'] = data['assignedReads']/data['assignedReads'].mean()
            data = pd.concat([data, d])
            data['exp'] = exp
            datas.append(data)

        self.data = pd.concat(datas).reset_index()

        self.widedata = self.data.pivot(index='#name', columns='exp', values='x')
        self.widedata.fillna(0, inplace=True)

        stat_data = self.data.pivot(index='#name', columns='exp', values='assignedReads')
        mean = stat_data.mean(axis=0, numeric_only=True)
        total = stat_data.sum(axis=0, numeric_only=True)
        
        stats = pd.DataFrame.from_dict({
            'mean': mean,
            'total': total
        })
        display(stats)