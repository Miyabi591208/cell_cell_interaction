import pandas as pd
import sys
import os

def correct_lr(data):
    '''
    Invert the RL to LR and R1R2 to r2>r1
    '''
    import pandas as pd
    def swap(a,b): return b,a
    data = data.to_dict('index')
    for k,v in data.items():
        if v['isReceptor_fst'] and v['isReceptor_scn']:
            v['isReceptor_fst'],v['isReceptor_scn'] = swap(v['isReceptor_fst'],v['isReceptor_scn'])
            v['Ligand'],v['Receptor'] = swap(v['Ligand'],v['Receptor'])
            v['Ligand.Cluster'],v['Receptor.Cluster'] = swap(v['Ligand.Cluster'],v['Receptor.Cluster'])
        elif v['isReceptor_fst'] and not v['isReceptor_scn']:
            v['isReceptor_fst'],v['isReceptor_scn'] = swap(v['isReceptor_fst'],v['isReceptor_scn'])
            v['Ligand'],v['Receptor'] = swap(v['Ligand'],v['Receptor'])
            v['Ligand.Cluster'],v['Receptor.Cluster'] = swap(v['Ligand.Cluster'],v['Receptor.Cluster'])
    res_df = pd.DataFrame.from_dict(data,orient='index')
    return (res_df)

def cpdb2df(data):
    data = data.fillna(0)
    df_data = {}
    df_data['Ligand'] = []
    df_data['Receptor'] = []
    df_data['Ligand.Cluster'] = []
    df_data['Receptor.Cluster'] = []
    df_data['isReceptor_fst'] = []
    df_data['isReceptor_scn'] = []
    df_data['MeanLR'] = []
    for i in range(data.shape[0]):
        pair = list(data['interacting_pair'])[i].split('_')
        for j in range(data.iloc[:,14:].shape[1]):
            c_pair = list(data.columns)[j+14].split('|')
            if float(data.iloc[i,j+14]) !=0.0:
                df_data['Ligand'].append(pair[0])
                df_data['Receptor'].append(pair[1])
                df_data['Ligand.Cluster'].append(c_pair[0])
                df_data['Receptor.Cluster'].append(c_pair[1])
                df_data['isReceptor_fst'].append(list(data['receptor_a'])[i])
                df_data['isReceptor_scn'].append(list(data['receptor_b'])[i])
                df_data['MeanLR'].append(data.iloc[i,j+14])
    data_final = pd.DataFrame.from_dict(df_data)
    return(data_final)

def cpdb2df_nocls(data):
    '''
   		When the cluster name is used on CPDB
    '''
    data = data.fillna(0)
    df_data = {}
    df_data['Ligand'] = []
    df_data['Receptor'] = []
    df_data['Ligand.Cluster'] = []
    df_data['Receptor.Cluster'] = []
    df_data['isReceptor_fst'] = []
    df_data['isReceptor_scn'] = []
    df_data['MeanLR'] = []
    for i in range(data.shape[0]):
        pair = list(data['interacting_pair'])[i].split('_')
        for j in range(data.iloc[:,12:].shape[1]):
            c_pair = list(data.columns)[j+12].split('|')
            if float(data.iloc[i,j+12]) !=0.0:
                df_data['Ligand'].append(pair[0])
                df_data['Receptor'].append(pair[1])
                df_data['Ligand.Cluster'].append(c_pair[0])
                df_data['Receptor.Cluster'].append(c_pair[1])
                df_data['isReceptor_fst'].append(list(data['receptor_a'])[i])
                df_data['isReceptor_scn'].append(list(data['receptor_b'])[i])
                df_data['MeanLR'].append(data.iloc[i,j+12])
    data_final = pd.DataFrame.from_dict(df_data)
    return(data_final)

# obtain argument
base_path = sys.argv[1]
path_dir = sys.argv[2]

# Check if file exists
if not os.path.exists(path_dir):
    raise FileNotFoundError(f"The file {path_dir} does not exist.")

# load data
s1 = pd.read_table(
    os.path.join(path_dir), 
    sep='\t'
)

s1_filtered = cpdb2df(s1)
s1_filtered = correct_lr(s1_filtered)

##### Adapting to new format
s1_filtered_final = pd.DataFrame()
s1_filtered_final['source'] = s1_filtered['Ligand.Cluster']
s1_filtered_final['target'] = s1_filtered['Receptor.Cluster']
s1_filtered_final['gene_A'] = s1_filtered['Ligand']
s1_filtered_final['gene_B'] = s1_filtered['Receptor']
s1_filtered_final['type_gene_A'] = 'Ligand'
s1_filtered_final['type_gene_B']= 'Receptor'
s1_filtered_final['MeanLR']=s1_filtered['MeanLR']

# save
s1_filtered_final.to_csv(
    os.path.join(base_path, "filtered_corrected.csv")
)
print(f"File saved Successfully")

