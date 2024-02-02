
import biomart                                                                                                                                                                                              

def get_ensembl_mappings():                                                                                                                                                                                 
    # Set up connection to server                                                                                                                                                                           
    server = biomart.BiomartServer('http://uswest.ensembl.org/biomart')                                                                                                                                     
    mart = server.datasets['mmusculus_gene_ensembl']                                                                                                                                                        

    # List the types of data we want                                                                                                                                                                        
    attributes = ['ensembl_transcript_id', 'mgi_symbol',                                                                                                                                                    
                  'ensembl_gene_id', 'ensembl_peptide_id', 'entrezgene_id']                                                                                                                                                  

    # Get the mapping between the attributes                                                                                                                                                                
    response = mart.search({'attributes': attributes})                                                                                                                                                      
    data = response.raw.data.decode('ascii')                                                                                                                                                                

    ensembl_to_genesymbol = {}                                                                                                                                                                              
    # Store the data in a dict                                                                                                                                                                              
    for line in data.splitlines():                                                                                                                                                                          
        line = line.split('\t')                                                                                                                                                                             
        # The entries are in the same order as in the `attributes` variable                                                                                                                                 
        transcript_id = line[0]                                                                                                                                                                             
        gene_symbol = line[1]                                                                                                                                                                               
        ensembl_gene = line[2]                                                                                                                                                                              
        ensembl_peptide = line[3]                                                                                                                                                                           

        # Some of these keys may be an empty string. If you want, you can
        # avoid having a '' key in your dict by ensuring the
        # transcript/gene/peptide ids have a nonzero length before
        # adding them to the dict
        ensembl_to_genesymbol[transcript_id] = gene_symbol
        ensembl_to_genesymbol[ensembl_gene] = gene_symbol
        ensembl_to_genesymbol[ensembl_peptide] = gene_symbol

    return ensembl_to_genesymbol
