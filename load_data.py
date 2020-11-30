import os
import gzip
from graphio import NodeSet, RelationshipSet, Container
from py2neo import Graph

from helper import GffReader


def parse_ncbigene(ncbigene_file):
    """
    Parser for gene_info.gz
    """

    # define a few things that are usually passed as parameter
    taxid = '9606'

    genes = NodeSet(['Gene'], merge_keys=['sid'])
    genesymbols = NodeSet(['GeneSymbol'], merge_keys=['sid', 'taxid'])
    genesymbol_synonym_genesymbol = RelationshipSet('SYNONYM', ['GeneSymbol'], ['GeneSymbol'],
                                                    ['sid', 'taxid'], ['sid', 'taxid'])
    gene_maps_genesymbol = RelationshipSet('MAPS', ['Gene'], ['GeneSymbol'], ['sid'], ['sid', 'taxid'])

    # check sets
    check_ids = set()
    check_ids_symbols = set()

    with gzip.open(ncbigene_file, 'rt') as f:

        header = next(f)

        # account for different formatting of header line (leading "#Format: " or not)
        if header.startswith('#Format:'):
            header_fields = tuple(header.split(':')[1].split('(')[0].rstrip().lstrip().split())
        elif header.startswith('#tax'):
            header_fields = tuple(
                header[1:].strip().split('\t')
            )
        else:
            raise AttributeError("File header was reformatted: {0}".format(header))

        for l in f:

            flds = l.rstrip().split('\t')

            this_taxid = flds[0]
            if this_taxid == taxid:

                # (Gene)
                entrez_gene_id = flds[1]
                if entrez_gene_id not in check_ids:
                    props = {'sid': entrez_gene_id, 'source': 'ncbigene'}
                    # update with all fields
                    props.update(
                        zip(header_fields, flds)
                    )

                    check_ids.add(entrez_gene_id)
                    genes.add_node(props)

                # (GeneSymbol) and (GeneSymbol)-[SYNONYM]-(GeneSymbol)
                primary_symbol = flds[2]
                synonym_symbols = flds[4].split('|')

                # add primary symbol node
                if primary_symbol not in check_ids_symbols and primary_symbol != '-':
                    check_ids_symbols.add(primary_symbol)
                    genesymbols.add_node({'sid': primary_symbol,
                                          'taxid': taxid})

                for synonym in synonym_symbols:
                    # GeneSymbol-[SYNONYM]-GeneSymbol
                    genesymbol_synonym_genesymbol.add_relationship({'sid': synonym, 'taxid': taxid},
                                                                   {'sid': primary_symbol,
                                                                    'taxid': taxid},
                                                                   {'source': 'ncbigene'})

                    if synonym not in check_ids_symbols and synonym != '-':
                        check_ids_symbols.add(synonym)
                        genesymbols.add_node({'sid': synonym,
                                              'status': 'synonym',
                                              'taxid': taxid})

                # (Gene)-[MAPS]-(GeneSymbol)
                # primary
                gene_maps_genesymbol.add_relationship({'sid': entrez_gene_id},
                                                      {'sid': primary_symbol, 'taxid': taxid},
                                                      {'source': 'ncbigene', 'status': 'primary'})
                # synonym
                for symbol in synonym_symbols:
                    gene_maps_genesymbol.add_relationship({'sid': entrez_gene_id},
                                                          {'sid': symbol, 'taxid': taxid},
                                                          {'source': 'ncbigene', 'status': 'synonym'})

    # add NodeSets/RelationshipSets to a container and return
    output = Container()
    output.add_all([genes, genesymbols, gene_maps_genesymbol, genesymbol_synonym_genesymbol])

    return output


def parse_ensembl(ensembl_gtf_file):
    """
    Parser for Homo_sapiens.GRCh38.101.chr_patch_hapl_scaff.gtf.gz
    """
    # define a few things that are usually passed to the function
    datasource_name = 'ensembl'
    taxid = '9606'

    # NodeSets
    genes = NodeSet(['Gene'], merge_keys=['sid'])
    transcripts = NodeSet(['Transcript'], merge_keys=['sid'])
    proteins = NodeSet(['Protein'], merge_keys=['sid'])

    # RelationshipSets
    gene_codes_transcript = RelationshipSet('CODES', ['Gene'], ['Transcript'], ['sid'], ['sid'])
    transcript_codes_protein = RelationshipSet('CODES', ['Transcript'], ['Protein'], ['sid'], ['sid'])

    annotation = GffReader(ensembl_gtf_file)

    check_gene_ids = set()
    check_transcript_ids = set()
    check_protein_ids = set()
    check_gene_transcript_rels = set()
    check_transcript_protein_rels = set()

    for r in annotation.records:

        # add gene node
        gene_id = r.attributes['gene_id']
        if gene_id not in check_gene_ids:
            props = {'sid': gene_id, 'name': r.attributes['gene_name'], 'taxid': taxid,
                     'source': datasource_name}

            genes.add_node(props)
            check_gene_ids.add(gene_id)

        # add transcript node
        if r.type == 'transcript':
            transcript_id = r.attributes['transcript_id']
            if transcript_id not in check_transcript_ids:
                props = {'sid': transcript_id, 'source': datasource_name, 'taxid': taxid}

                transcripts.add_node(props)
                check_transcript_ids.add(transcript_id)

        # add protein node
        if r.type == 'CDS':
            protein_id = r.attributes['protein_id']
            if protein_id not in check_protein_ids:
                props = {'sid': protein_id, 'source': datasource_name, 'taxid': taxid}

                proteins.add_node(props)
                check_protein_ids.add(protein_id)

        # Gene-CODES-Transcript
        if r.type == 'transcript':
            transcript_id = r.attributes['transcript_id']
            gene_id = r.attributes['gene_id']

            # add gene-transcript rel
            if gene_id + transcript_id not in check_gene_transcript_rels:
                gene_codes_transcript.add_relationship({'sid': gene_id}, {'sid': transcript_id},
                                                       {'source': datasource_name})
                check_gene_transcript_rels.add(gene_id + transcript_id)

        # Transcript-CODES-Protein
        if r.type == 'CDS':
            protein_id = r.attributes['protein_id']
            transcript_id = r.attributes['transcript_id']

            # add transcript-protein rel
            if transcript_id + protein_id not in check_transcript_protein_rels:
                transcript_codes_protein.add_relationship({'sid': transcript_id}, {'sid': protein_id},
                                                          {'source': datasource_name})
                check_transcript_protein_rels.add(transcript_id + protein_id)

    output = Container()
    output.add_all([genes, transcripts, proteins, gene_codes_transcript, transcript_codes_protein])
    return output


def parse_ensembl_ncbi_mappings(ensembl_tsv_file):
    """
    Second parser for Homo_sapiens.GRCh38.101.entrez.tsv.gz
    """
    # define a few things that are usually passed to the function
    datasource_name = 'ensembl'
    taxid = '9606'

    gene_maps_gene = RelationshipSet('MAPS', ['Gene'], ['Gene'], ['sid'], ['sid'])

    check_rels = set()

    with gzip.open(ensembl_tsv_file, 'rt') as f:
        lines = f.readlines()
        for l in lines[1:]:
            flds = l.strip().split()

            ensembl_gene_id = flds[0]
            ncbi_gene_id = flds[3]

            if frozenset([ensembl_gene_id, ncbi_gene_id]) not in check_rels:
                gene_maps_gene.add_relationship({'sid': ensembl_gene_id}, {'sid': ncbi_gene_id},
                                                     {'source': datasource_name})
                check_rels.add(frozenset([ensembl_gene_id, ncbi_gene_id]))

    output = Container()
    output.add_all([gene_maps_gene])
    return output


# graph
graph = Graph(host='localhost', user='neo4j', password='test', name='graphiopy2neo')

# get data files
this_path = os.path.dirname(os.path.abspath(__file__))
ncbigene_file = os.path.join(this_path, 'Homo_sapiens.gene_info.gz')
ensembl_gtf_file = os.path.join(this_path, 'Homo_sapiens.GRCh38.101.chr_patch_hapl_scaff.gtf.gz')
ensembl_tsv_file = os.path.join(this_path, 'Homo_sapiens.GRCh38.101.entrez.tsv.gz')

# run parser
ncbi_gene_data = parse_ncbigene(ncbigene_file)
ensembl_data = parse_ensembl(ensembl_gtf_file)
ensembl_ncbi_mapping_data = parse_ensembl_ncbi_mappings(ensembl_tsv_file)

# load data to Neo4j
## start with nodes
for datacontainer in [ncbi_gene_data, ensembl_data, ensembl_ncbi_mapping_data]:
    for nodeset in datacontainer.nodesets:
        nodeset.create_index(graph)
        nodeset.create(graph)

## then relationships
for datacontainer in [ncbi_gene_data, ensembl_data, ensembl_ncbi_mapping_data]:
    for relationshipset in datacontainer.relationshipsets:
        relationshipset.create_index(graph)
        relationshipset.create(graph)
