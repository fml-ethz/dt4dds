import numpy as np
import pandas as pd
import pathlib
import ruamel.yaml
import Bio.SeqIO
import Bio.SeqUtils
import re

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
        data = ruamel.yaml.safe_load(open(file))
        self.raw_data = {}
        self.data = {}
        self.exp = self.file.parent.name

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

        
        # OVERALL ERROR RATES
        df = pd.DataFrame({
            'type': ('substitutions', 'deletions', 'insertions', 'delevents'),
            'rate': [getattr(self, f"r_{errortype}") for errortype in ('substitutions', 'deletions', 'insertions', 'delevents')]
        })
        df['read'] = str(self.read)
        df['exp'] = self.exp
        self.data[f"overall_error_rates"] = df

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
            dfs.append(df)
        self.data[f"overall_error_rates_by_position"] = pd.concat(dfs).reset_index(drop=True)

        # COMBINED ERROR RATES BY REFPOSITION
        dfs = []
        for errortype in ('substitutions', 'deletions', 'insertions'):
            df = self.data[f"{errortype}_by_refposition"].copy()
            df['type'] = errortype
            dfs.append(df)
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
                dfs.append(df)
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
                dfs.append(df)
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
            dfs.append(df)
        self.data[f"error_frequency_by_length"] = pd.concat(dfs).reset_index(drop=True)

        # ERROR FREQUENCY PER READ
        dfs = []
        for errortype in ('substitutions', 'delevents', 'insertions'):
            df = self.raw_data[f"p_{errortype}_by_read"].reset_index()
            df = df.rename(columns={'index': 'frequency'})
            df['type'] = errortype
            df['read'] = str(self.read)
            df['exp'] = self.exp
            dfs.append(df)
        self.data[f"error_frequency_by_read"] = pd.concat(dfs).reset_index(drop=True)






class ErrorAnalysis():
    
    def __init__(self, folder, skip=[(), ()], paired=False, group="matches", suffix="stats", r1_file="R1.fq.gz", r2_file="R2.fq.gz"):

        # build the file path to the analysis files
        paired_suffix = "paired." if paired else ""
        r1_path = f"{r1_file}.{paired_suffix}{group}.{suffix}"
        r2_path = f"{r2_file}.{paired_suffix}{group}.{suffix}"
        self.folder = pathlib.Path(folder)

        # read analysis files
        self.r1 = ErrorFile(self.folder / r1_path, read=1, skip=skip[0])
        self.r2 = ErrorFile(self.folder / r2_path, read=2, skip=skip[1])
        reads = [self.r1, self.r2]

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






MIN_HOMOPOLYMER_LENGTH = 2

class SeqErrorAnalysis():

    def __init__(self, folder, paired=True):

        self.folder = pathlib.Path(folder)
        self.exp = self.folder.name
        self.paired = paired
        self.ref_dict = read_reference(self.folder / "design_files.fasta")

        # HOMOPOLYMER DATA
        self.homopolymer_dict = {}
        max_homopolymer_dict = {}
        for refid, seq in self.ref_dict.items():
            self.homopolymer_dict[refid] = []
            for base in ('A', 'C', 'G', 'T'):
                for match in re.compile(fr"{base*MIN_HOMOPOLYMER_LENGTH}+").finditer(str(seq.seq)):
                    self.homopolymer_dict[refid].append((base, match.end()-match.start(), match.start(), match.end()))
            max_homopolymer_dict[refid] = max(self.homopolymer_dict[refid], key=lambda x: x[1])[1]
        homopolymer_data = pd.DataFrame.from_dict(max_homopolymer_dict, orient="index", columns=['max_HP']).reset_index().rename(columns={'index': 'refid'})


        # ABUNDANCE DATA
        abundance_data = pd.read_csv(self.folder / "scafstats.txt", sep="\t", converters={'#name': str})
        abundance_data.rename(columns={"#name": "refid"}, inplace=True)
        abundance_data['rel_abundance'] = abundance_data['assignedReads']/abundance_data['assignedReads'].mean()
        abundance_data = abundance_data[['refid', 'rel_abundance']]

        # FREE ENERGY DATA
        energy_data = pd.read_csv(self.folder / "secondary_energy.csv", header=0, names=["refid", "energy_fw", "energy_rv"], converters={'refid': str})
        energy_data = energy_data[['refid', 'energy_fw', 'energy_rv']]

        dfs = []
        self.mean_rates = {}
        n_reads = 2 if paired else 1
        for read in range(1, n_reads+1):
            self.mean_rates[read] = {}
            df = pd.read_csv(self.folder / f"R{read}.fq.gz.overview_by_sequence.csv", converters={'refid': str})
            df['read'] = str(read)
            df['exp'] = self.exp
            df['GC'] = [Bio.SeqUtils.gc_fraction(self.ref_dict[refid].seq) for refid in df['refid']]

            for e in ('substitutions', 'deletions', 'insertions'):
                df[f'r_{e}'] = df[f'n_{e}']/df['n_bases']

            for e in ('substitutions', 'deletions', 'insertions'):
                self.mean_rates[read][e] = df[f'n_{e}'].sum()/df['n_bases'].sum()
                df[f'rel_{e}'] = df[f'r_{e}']/self.mean_rates[read][e]
                df.drop([f'n_{e}'], axis=1, inplace=True)
            df.drop([f'n_bases'], axis=1, inplace=True)
            
            df = df.merge(homopolymer_data, on="refid")
            df = df.merge(abundance_data, on="refid")
            df = df.merge(energy_data, on="refid")
            dfs.append(df)
        self.overview_data = pd.concat(dfs).reset_index(drop=True)

        self.position_data = {}
        for read in range(1, n_reads+1):
            self.position_data[read] = {}
            for errortype in ('substitutions', 'deletions', 'insertions'):
                df = pd.read_csv(self.folder / f"R{read}.fq.gz.{errortype}_by_sequence.csv", converters={'refid': str}).set_index("refid")
                self.position_data[read][errortype] = df



    def prepare_seqanalysis(self):

        overview = self.overview_data.drop(self.overview_data[self.overview_data['read'] != "1"].index).set_index("refid")
        data = {
            'exp': [],
            'refid': [],
            'length': [],
            'n_reads': [],
            'type': [],
            'rate': [],
            'rel_rate': []
        }

        for refid in overview.index.unique():

            n_reads = overview.loc[refid, "n_reads"]
            substitution_window = self.position_data[1]['substitutions'].loc[refid]
            deletion_window = self.position_data[1]['deletions'].loc[refid]
            insertion_window = self.position_data[1]['insertions'].loc[refid]

            hp_list = self.homopolymer_dict.get(refid, [])
            to_skip = set()

            for base, length, start, end in hp_list:
                to_skip.update(range(start, end))
                rates = {}
                rates['substitutions'] = np.average(substitution_window[start:end])
                rates['deletions'] = np.average(deletion_window[start:end])
                rates['insertions'] = np.average(insertion_window[start:end])

                for type, r in rates.items():
                    data['exp'].append(self.exp)
                    data['refid'].append(refid)
                    data['length'].append(length)
                    data['n_reads'].append(n_reads)
                    data['type'].append(type)
                    data['rate'].append(r)
                    data['rel_rate'].append(r/self.mean_rates[1][type])

            for errortype, errorwindow in (('substitutions', substitution_window), ('deletions', deletion_window), ('insertions', insertion_window)):
                to_include = set(range(len(errorwindow))).difference(to_skip)
                rate = np.average(errorwindow[list(to_include)])

                data['exp'].append(self.exp)
                data['refid'].append(refid)
                data['length'].append(1)
                data['n_reads'].append(n_reads)
                data['type'].append(errortype)
                data['rate'].append(rate)
                data['rel_rate'].append(rate/self.mean_rates[1][errortype])

        df = pd.DataFrame.from_dict(data)
        df.to_csv(self.folder / "seqanalysis.csv", index=False)




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
            refset = set(read_reference(pathlib.Path(file).parent / "design_files.fasta").keys())
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