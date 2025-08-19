from Bio import Entrez
from Bio import SeqIO
import pandas as pd


def false_genes(email_address, organism = None):
    Entrez.email = email_address
    if organism is None:
        search = 'Record[All Fields] AND support[All Fields] AND submission[All Fields] AND GeneRIFs[All Fields] AND a gene[All Fields]) AND Bacteria[porgn] AND alive[prop]'
    else:
        search = f'Record[All Fields] AND support[All Fields] AND submission[All Fields] AND GeneRIFs[All Fields] AND a gene[All Fields]) AND {organism}[porgn] AND alive[prop]'
    record = Entrez.read(Entrez.esearch(db="gene",term=search, retmax = 10000)) #porgn - primary organism
    id_list = record["IdList"]
    return(id_list)

def retrieve_pub_from_gene(gene_nm, email_address, rm_false = [], organism = None):
    """
    Function retrieves all publications associated with a gene search, the user can input 4 arguments:
     1) gene_nm: could be a gene name (eg. PhoR) or geneID or its full name (eg. sensor histidine kinase PhoR)
     2) email: your personal email
     3) rm_false: can choose to add in list of genes that are false placeholders, default is an empty list
     4) organism: can specify to search for genes within species, genus, phylum, etc, default is to search for the gene within all bacteria

    Returns a list and a dataframe stored within one object, first is a list of publications associated with gene of interest, second is dataframe with one column of geneIDs, and a second column showing list of publciations associated with each gene
    """
    Entrez.email = email_address
    if organism is None:
        search = 'Bacteria[porgn]'+ gene_nm + ' AND alive[prop]'
    else:
        search = f'{organism}[porgn]' + gene_nm + ' AND alive[prop]'
    record = Entrez.read(Entrez.esearch(db="gene",term=search, retmax = 10000)) #porgn - primary organism
    #record = Entrez.read(handle)
    id_list = record["IdList"]
    if len(rm_false) > 0:
        id_list = [item for item in id_list if item not in rm_false]
    pub_record = Entrez.read(Entrez.elink(db = "pubmed", dbfrom="gene", id=id_list))
    pub_list = []
    pub_per_gene = []
    for i in range(len(pub_record)):
        if len(pub_record[i]['LinkSetDb']) > 0:
            link_index = next((i for i, item in enumerate(pub_record[i]["LinkSetDb"]) if item.get('LinkName') == 'gene_pubmed'), -1)
            pub_per_gene.append({'geneID': pub_record[i]['IdList'][0], 'num_publications': len(pub_record[i]['LinkSetDb'][link_index]['Link'])})
 
            if link_index != -1:
                for link in pub_record[i]["LinkSetDb"][link_index]["Link"]:
                    pub_list.append(link["Id"])
    
    pub_per_gene_df = pd.DataFrame(pub_per_gene).sort_values(by = ['num_publications'], ascending=False)
    return([pub_list, pub_per_gene_df])


def retrieve_pub_from_term(search_key, email_address, search_field = ''):
    Entrez.email = email_address
    handle = Entrez.esearch(db="pubmed",term= search_key + f'{search_field}', retmax = 10000) 
    record = Entrez.read(handle)
    id_list = record["IdList"]
    return(id_list)

def gene_pub_search(gene_nm, search_key, email_address, organism = None, search_field = ''):
    gene_search = retrieve_pub_from_gene(gene_nm, email_address, organism)
    pub_search = retrieve_pub_from_term(search_key, email_address, search_field)
    intersect_list = list(set(gene_search[0]) & set(pub_search))
    return(intersect_list)

# can specify filters by organism 


def retrieve_seq_IDs(search_term, search_db, organism, max_len, email_address):
    Entrez.email = email_address
    if search_db == 'gene':
        IDs = Entrez.read(Entrez.esearch(db=search_db, term= f'{search_term}' + " AND " + f'{organism}[porgn]', retmax = 10000))["IdList"]    
    elif search_db in ('nucleotide','protein'):
        IDs = Entrez.read(Entrez.esearch(db=search_db, term= f'{search_term}' + " AND " + f'{organism}[porgn] ' + "AND (\"1\"[SLEN]:\"" + f'{max_len}' + "\"[SLEN])", retmax = 10000))["IdList"]
    else:
        print("Please enter 'gene', 'nucleotide', or 'protein' for the search_db parameter")
        return

    num_retstart = 10000
    num_IDs = len(IDs)
    if num_IDs == 10000:

        if search_db == 'gene':
            add_IDs = Entrez.read(Entrez.esearch(db=search_db, term= f'{search_term}' + " AND " + f'{organism}[porgn]', retmax = 10000, retstart = num_retstart))["IdList"]
        else:
            add_IDs = Entrez.read(Entrez.esearch(db=search_db, term= f'{search_term}' + " AND " + f'{organism}[porgn] '+"AND (\"1\"[SLEN]:\"" + f'{max_len}' + "\"[SLEN])", retmax = 10000, retstart = num_retstart))["IdList"]
        IDs = IDs + add_IDs
        num_retstart = num_retstart + 10000
        num_IDs = len(add_IDs)
    return(IDs)
    

# retrieve fasta sequences by search term
def retrieve_fasta_seqs(search_term, search_db, organism, max_len, email_address, filename):
    """
    Function retrieves all nucleotide sequences associated with a search term, the user can input 5 arguments:
     1) search term: could be a gene name (eg. PhoR) or its full name (eg. sensor histidine kinase PhoR)
     2) search_db: database to search against, can be gene, nucleotide, or protein 
     2) organism: only retrieve sequences belonging to this organism (can be species, family, order, etc), can be a list of organisms or just the string name of one organism
     3) max_len: maximum length of sequence, setting this parameter can get rid of whole genome sequences
     4) email: your personal email
     5) filename: name of the output file

    A file with all fasta nucleotide sequences associated with search term is generated 
    """
    Entrez.email = email_address
    if isinstance(organism, list):
        IDs = []
        for species in organism:
            IDs = IDs + retrieve_seq_IDs(search_term = search_term, search_db = search_db, organism = species, max_len = max_len, email_address = email_address)
    else:
        IDs = retrieve_seq_IDs(search_term = search_term, search_db = search_db, organism = organism, max_len = max_len, email_address = email_address)
    
    with open(filename, 'w') as f:
        print(f'{search_db} ' + "fasta sequences for " + f'{organism}\n',file = f)
        
    with open(filename, 'a') as f:
        print(Entrez.efetch(db=search_db, id=IDs, rettype="fasta", retmode="text").read(), file = f)


def retrieve_associated_tax(search_term, search_db, organism, max_len, email_address):
    Entrez.email = email_address
    IDs = retrieve_seq_IDs(search_term, search_db, organism, max_len, email_address)
    print('getting IDs')

    handle = Entrez.efetch(db=search_db, id=IDs, retmode = "xml")
    records = Entrez.read(handle)
    print(len(records))
    handle.close()
    tax_list = []
    for record in records:
        if search_db == 'gene':
            if "Gene-track_status" in record['Entrezgene_track-info']['Gene-track'] and record['Entrezgene_track-info']['Gene-track']['Gene-track_status'].attributes['value'] == 'discontinued':
                continue
            tax_list.append(record['Entrezgene_source']['BioSource']['BioSource_org']['Org-ref']['Org-ref_taxname'])
        else:
            tax_list.append(record['GBSeq_organism'])

    final_taxonomy = list(set(tax_list))
    return(final_taxonomy)

