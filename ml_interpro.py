#!/usr/bin/env python
print('interpro IPR001578')
#import topoly
import requests, sys, json
import csv
import urllib
import re
import ast
import pandas as pd
import operator
import itertools
import math
import copy

WEBSITE_API = "https://rest.uniprot.org/beta"

def get_url(url, **kwargs):
    response = requests.get(url, **kwargs)
    
    if not response.ok:
        print(response.text)
        response.raise_for_status()
        sys.exit()
        
    return response

def uniprot_structures(ID):
    result = {}
    data = get_url(f"{WEBSITE_API}/uniprotkb/{ID}?fields=xref_alphafolddb, xref_pdb, sequence")
    data = data.text
    pdb = re.compile(r"\"PDB\",\"id\":\"(\w+)\"")
    pdbs = ' '.join(re.findall(pdb, data))
    result["PDB"] = pdbs
    af = re.compile(r"\"AlphaFoldDB\",\"id\":\"(\w+)\"")
    afs = ' '.join(re.findall(af, data))
    result["AF"] = afs
    seq = re.compile(r"\"sequence\":{\"value\":\"(\w+)\"")
    seqs = ''.join(re.findall(seq, data))
    result["sequence"] = seqs
    return result

def get_cif(ID):
    text = urllib.request.urlopen('https://alphafold.ebi.ac.uk/files/AF-'+ID+'-F1-model_v3.cif').read().decode('utf-8')
    pattern = re.compile(r"data_(AF-\w+-\w+)")
    matches = re.findall(pattern, text)
    return text, matches[0]

def get_cif_ID(ID):
    text = urllib.request.urlopen('https://alphafold.ebi.ac.uk/files/AF-'+ID+'-F1-model_v3.cif').read().decode('utf-8')
    pattern = re.compile(r"data_(AF-\w+-\w+)")
    matches = re.findall(pattern, text)
    open(ID+".cif", "w+").write(text)

def get_pLDDT(cif):
    pattern = re.compile(r"_ma_qa_metric_global\.metric_value\s(\d+\.\d+)")
    match = float(re.findall(pattern, cif)[0])
    return match

def get_IDs(clstr_file):
    result = []
    text = open(clstr_file).read()
    clusters = text.split(">Cluster")[1:]
    clusters = [[res.split('>') for res in cluster.split('\n')[1:-1]] for cluster in clusters]
    for i in range(len(clusters)):
        dct = {}
        dct["nr"] = i
        dct["rest"] = []
        p_length = re.compile(r"\d+\s+(\d+)aa,")
        star_seq = re.compile(r"(\w+)[A-Za-z0-9|,\-().;/:+_]+\s\*")
        for struct in clusters[i]:
            if re.findall(star_seq, struct[1]):
                ID = re.findall(star_seq, struct[1])[0]
                dct["representative"] = (ID, int(re.findall(p_length, struct[0])[0]))
            else:
                pattern = re.compile(r"(\w+)[A-Za-z0-9|,\-().;/:+_]+\s")
                ID = re.findall(pattern, struct[1])[0]
                try:
                    if ID != dct["representative"]:
                        dct["rest"].append((ID,int(re.findall(p_length, struct[0])[0])))
                except:
                    dct["rest"].append((ID, int(re.findall(p_length, struct[0])[0])))
        result.append(dct)
    return result

def count_knot(cif_file_content):
    open("cif.cif", "w+").write(cif_file_content)
    #return topoly.homfly("cif.cif", tries=400, max_cross=30, translate=True)
    return topoly.alexander("cif.cif", tries=500, max_cross=30, translate=True)


def knot_matrix(cif_file_name, seq_length):
    if seq_length >= 5000:
        params = ALEX_PARAMS[5000]
    elif seq_length >= 2500:
        params = ALEX_PARAMS[2500]
    elif seq_length >= 1500:
        params = ALEX_PARAMS[1500]
    elif seq_length >= 1000:
        params = ALEX_PARAMS[1000]
    else:
        params = ALEX_PARAMS[250]

    return find_cores(topoly.alexander(cif_file_name, matrix=True, max_cross=30, 
            matrix_density = params[0], tries=25, matrix_calc_cutoff = 0.2, translate=True)) 

def save_clusters(clusters_filename, clstr_file):
    with open(clusters_filename, 'w+') as file:
        clusters = get_IDs(clstr_file)
        file.write(str(clusters))

#save_clusters("IPR013128_clusters.txt", "IPR013128_peptidase_c1a_representants095.fasta.clstr")

def choose_knot(topoly_output):
    print("Topoly output:",topoly_output)
    if "0_1" in topoly_output.keys(): # if trivial knot is in the set of keys
        if topoly_output["0_1"] >= 0.5: # if probability of unknot is greater or equal to 50% we decide that it is actually unknotted
            knot = "0_1"
            prob = topoly_output["0_1"]
        else: # if probability of unknot is lower than 50% we want to check if there is one dominant knot with at least probability 30%
            topoly_out = copy.deepcopy(topoly_output)
            topoly_output.pop('0_1', None)
            knot, prob =  max(topoly_output.items(), key=operator.itemgetter(1)) # we take the knot with the highest probability
            if prob < 0.35: # if the probability is lower than 30% we cannot assign this protein knotted topology
                knot, prob = "0_1", topoly_out["0_1"] # in this case we assign unknot even if the probability of unknot is lower than 50%
    else: # if unknot is not in topoly output we want to assign knotted topology
        knot, prob =  max(topoly_output.items(), key=operator.itemgetter(1))
        if prob < 0.35: # but if the highest probability is lower than 30% we cannot assign this knot so we assign None
            knot, prob = None, None # might look strange in results but this case is unlikely to happen
    print("Chosen topology and probability", knot, prob)
    return (knot, prob)

def write_info(out_filename, clusters_file):
    with open(out_filename,'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['Uniprot ID', 'sequence length', 'No. cluster', 'sequence',
                             'PDB IDs', 'representative structure', 'AlphaFold ID', 'pLDDT', 
                             'knot', 'knot probability', 'topoly output'])
        clusters = ast.literal_eval( open(clusters_file).read())
        for cluster in clusters:
            #done = False
            for ID, length in [cluster["representative"]]+cluster["rest"]:
                try:
                    cif, af_id = get_cif(ID)
                    plddt = get_pLDDT(cif)
                    if plddt > 70:
                        topology = count_knot(cif)
                        knot, prob = choose_knot(topology)
                        uniprot_data = uniprot_structures(ID)
                        tsv_writer.writerow([ID, length, cluster["nr"], 
                                        uniprot_data['sequence'], uniprot_data["PDB"], 
                                        ID==cluster["representative"][0], af_id, 
                                        plddt, knot, prob, topology])
                except urllib.error.HTTPError:
                    uniprot_data = uniprot_structures(ID)
                    tsv_writer.writerow([ID, length, cluster["nr"], 
                                        uniprot_data['sequence'], uniprot_data["PDB"], 
                                        ID==cluster["representative"][0], "", 
                                        "", "", "", ""])



ALEX_PARAMS = {
# length:density
250 : [3, 5],
1000: [7, 12],
1500: [7, 20],
2500: [9, 20],
5000: [9, 20]
}  

def save_structures(clusters_file):
    clusters = ast.literal_eval( open(clusters_file).read())
    for cluster in clusters:
        for ID, length in [cluster["representative"]]+cluster["rest"]:
            try:
                get_cif_ID(ID)
            except urllib.error.HTTPError:
                pass

def count_cores(PFAM_ID):
    # jeszcze trzeba sprawdzić czy tam w ogóle jest węzeł! bo jest źle zrobione
    # 0_1 > 0.5 nie powinno liczyć macierzy
    with open(PFAM_ID+'output_info.csv') as info_file:
        with open(PFAM_ID+'_cores_from_matrix_calculation.csv', 'wt') as out_file:
            tsv_writer = csv.writer(out_file, delimiter='\t')
            tsv_writer.writerow(['UNIPROT ID', 'pLDDT', 'knot', 'knot probalility', 'krot core'])
            df = pd.read_csv(info_file, sep='\t', header = 0)
            nrows = df.shape[0]
            for i in range(nrows):
                if df.loc[i, 'count matrix']:
                    tsv_writer.writerow([','.join(df.loc[i,['Uniprot ID','pLDDT', 'knot', 'knot probability']].to_string(header=False, index=False).split('\n')), 
                                        find_cores(knot_matrix('cif.cif', df.loc[i, ['sequence length']]))])





def find_cores(D, freq=0.42, size=2, dens=1):
    topologies = set([])

    cores = {} # {'3_1': (20,54), '4_1': ...}
    # remove empty elements of the dictionary
    D = {k: v for k, v in D.items() if v}
    for pos, topolprob in D.items():
        if topolprob:
            topol, prob = sorted(list(topolprob.items()), key=lambda x: x[1])[-1]
            if prob >= freq and topol!='0_1':
                topologies.add(topol)

                if topol in cores.keys() and (cores[topol][1]-cores[topol][0]) <= (pos[1]-pos[0]):
                    continue
                new_core = True
                for i, j in itertools.product(range(0,size*dens+1,dens), range(0,size*dens+1,dens)):
                    new_pos = (pos[0]-i,pos[1]+j)
                    if not (new_pos in D):
                        new_core = False
                        continue
                    topol2, prob2 = sorted(list(D[new_pos].items()), key=lambda x: x[1])[-1]
                    if topol2 != topol or prob2 < freq:
                        new_core = False
                        continue
                #if new_core: cores[topol] = [pos[1]-pos[0]+1,pos]
                if new_core: cores[topol] = pos
    return cores

write_info("IPR013128_peptidase_c1a_AF_info_new.csv", "IPR013128_clusters.txt")




