import os
os.chdir("Documents\\soil") # IMPORTANT: replace "Documents\\soil" with the path that leads to directory where the ncbi_scrape.py script is stored

from ncbi_scrape import *

# replace with your own email
user_email = 'angezou@gmail.com'

# First need to get list of placeholder genes (ie. genes with description "Record to support submission of GeneRIFs for a gene not in Gene") that have many publications attributed to them 
# The following command gets all placeholder genes in bacteria, can also get placeholder genes specific to a species, class, order, etc, just add organism = [organism of interest] (eg. 'Pseudomonas)
genes_to_rm = false_genes(user_email)

# The following function retrieve_pub_from_gene retrieves all publications associated with a gene search, the user can input 4 arguments:
# 1) gene_nm: could be a gene name (eg. PhoR) or geneID or its full name (eg. sensor histidine kinase PhoR)
# 2) email: your personal email
# 3) rm_false: can choose to add in list of genes that are false placeholders, default is an empty list
# 4) organism: can specify to search for genes within species, genus, phylum, etc, default is to search for the gene within all bacteria

# the following commands searches for the 'PhoR' gene in all of bacteria after removing the placeholder genes, command should take 5-10 min to run
phor_all = retrieve_pub_from_gene(gene_nm = 'PhoR', email_address = user_email, rm_false = genes_to_rm)

# phor_all stores two things, the first one is the list of publications associated with PhoR bacterial genes
phor_all[0]

# the second thing is a dataframe showing all the PhoR genes as listed by there geneID, and the # of publications associated with each gene 
phor_all[1].head(20)

# can also retrieve publicatiosn for all PhoR genes SPECIFIC to pseudomonas 
phor_pseudo = retrieve_pub_from_gene('PhoR', user_email, genes_to_rm, 'Pseudomonas')

# The following function retrieves all publications according to search term, can also limit search within a specific term (eg. Title, or abstract) by entering a 3rd argument eg. search_term = 'Title'
phor_pheno = retrieve_pub_from_term(search_key='phosphate solubilization bacteria', email_address=user_email)

#intersect the results of retrieve_pub_from_gene and retrieve_pub_from_term to get all publications with geneID AND phenotype
list(set(phor_all[0]) & set(phor_pheno))
#outputs 2 publications ['31405015', '32414048']

#Can perform all the previous steps to get publications with both gene association and phenotype data
phor_cross_ref = gene_pub_search(gene_nm = 'phosphate regulon sensor histidine kinase', search_key='phosphate solubilizing bacteria', user_email)
list(set(phor_all[0]) & set(phor_pheno))

#retriee all nucleotide sequences for nifH found in Bradyrhizobium elkanii
retrieve_fasta_seqs(search_term="nifh nitrogenase iron protein", 
                   search_db="nucleotide",
                   organism="Bradyrhizobium elkanii", 
                   max_len=1500,
                   email_address="angezou@gmail.com",
                   filename="nifh_Brady_output.txt")

#retrieve all nucleotide sequences for nifH found in Bradyrhizobium elkanii
retrieve_fasta_seqs(search_term="nifh nitrogenase iron protein", 
                   search_db="protein",
                   organism="Klebsiella", 
                   max_len=500, 
                   email_address="angezou@gmail.com",
                   filename="nifh_Klebsiella_output.txt")


#retrieve all bacterial species that contain nifh
nifh = retrieve_associated_tax(search_term="nifh nitrogenase iron protein", 
                   search_db="protein", #switch to "gene" if downloading nucleotide sequences
                   organism="Klebsiella variicola", 
                   max_len=800, 
                   email_address="angezou@gmail.com")

#retrieve all bacterial species that contain nifd
nifd = retrieve_associated_tax(search_term="nitrogenase molybdenum-iron protein alpha chain", 
                   search_db="protein", #switch to "gene" if downloading nucleotide sequences
                   organism="Klebsiella variicola", 
                   max_len=800, 
                   email_address="angezou@gmail.com")

#retrieve all bacterial species that contain nifk
nifk = retrieve_associated_tax(search_term="nitrogenase molybdenum-iron protein subunit beta", 
                   search_db="protein", #switch to "gene" if downloading nucleotide sequences
                   organism="Klebsiella variicola", 
                   max_len=800, 
                   email_address="angezou@gmail.com")

#get list of bacterial species that contain nifh, nifd, nifk (making them relevant nitrogen fixers)
nitrogen_fixers = list(set(nifh).intersection(nifd,nifk))

#retrieve fasta sequences for nifh given the list of nitrogen fixing organisms
retrieve_fasta_seqs(search_term="nifh nitrogenase iron protein", 
                   search_db="protein",
                   organism = nitrogen_fixers, 
                   max_len=800, 
                   email_address="angezou@gmail.com",
                   filename="nit_fix_output_protein.txt")

