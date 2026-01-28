from scripts.data_structure.wiki_data_structure import (
    Xref, Graphics, DataNode, DataNodeType, HAlign, VAlign,
    BorderStyle, ShapeType, Property, Comment, CitationRef, Citation,
    Annotation, AnnotationType, AnnotationRef
)
from scripts.parsing_functions import parsing_utils
from scripts.object2gmpl.gpml_writer import GPMLWriter
from scripts.utils.HTML_cleaner import clean_text_label
from scripts.utils import standard_graphics
from scripts.utils.property_parser import (
    create_dynamic_properties, handle_unique_id, handle_synonyms
)
import uuid
import re
from pathlib import Path


# Organism mapping cache
_ORGANISM_MAPPING = None

def load_organism_mapping():
    """Load organism mapping from org_id_mapping_v2.tsv."""
    global _ORGANISM_MAPPING

    if _ORGANISM_MAPPING is not None:
        return _ORGANISM_MAPPING

    _ORGANISM_MAPPING = {}

    # Load mapping file created at pipeline start
    mapping_file = Path(__file__).parent.parent.parent / 'org_id_mapping_v2.tsv'
    if mapping_file.exists():
        with open(mapping_file, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('org_id'):  # Skip header
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    organism_id = parts[0]
                    latin_name = parts[1]
                    if latin_name:
                        _ORGANISM_MAPPING[organism_id] = latin_name

    return _ORGANISM_MAPPING


def convert_to_latin_name(organism_id):
    """
    Convert ORG-ID or TAX-ID to Latin name.

    Args:
        organism_id: String like 'ORG-5993' or 'TAX-3702'

    Returns:
        str: Latin name if found, otherwise original organism_id
    """
    if not organism_id:
        return organism_id

    mapping = load_organism_mapping()
    organism_id = str(organism_id).strip()

    # Direct lookup (handles both ORG-XXXX and TAX-XXXX)
    if organism_id in mapping:
        return mapping[organism_id]

    # Return original if not found
    return organism_id


def create_species_annotation(record):
    """
    Create Annotation and AnnotationRef for species/taxonomy.

    Args:
        record: Gene/Protein record dictionary

    Returns:
        tuple: (list of Annotation objects, list of AnnotationRef objects)
    """
    annotations = []
    annotation_refs = []

    if 'SPECIES' in record:
        species_list = record['SPECIES']
        if not isinstance(species_list, list):
            species_list = [species_list]

        for species in species_list:
            species_id = str(species).strip()
            latin_name = convert_to_latin_name(species_id)
            
            # Create a unique ID for the annotation
            annotation_id = f"taxonomy_{species_id}"
            # Sanitize ID
            annotation_id = re.sub(r'[^a-zA-Z0-9_]', '_', annotation_id)
            
            annotation = Annotation(
                elementId=annotation_id,
                value=latin_name,
                type=AnnotationType.TAXONOMY,
                xref=Xref(identifier=species_id, dataSource="Taxonomy")
            )
            annotations.append(annotation)
            annotation_refs.append(AnnotationRef(elementRef=annotation_id))

    return annotations, annotation_refs


def parse_gene_dblinks(dblinks):
    """
    Parse DBLINKS to extract database references for genes

    Args:
        dblinks: List of database link strings

    Returns:
        dict: Dictionary with database names as keys and IDs as values
    """
    db_refs = {}

    for dblink in dblinks:
        if isinstance(dblink, str):
            try:
                match = re.match(r'\((\w+)\s+"([^"]+)"', dblink)
                if match:
                    db_name = match.group(1)
                    db_id = match.group(2)
                    db_refs[db_name] = db_id
            except Exception as e:
                print(f"Warning: Could not parse dblink: {dblink}")

    return db_refs


def get_primary_gene_xref(db_refs, record=None):
    """
    Get the primary external reference for a gene

    Args:
        db_refs: Dictionary of database references from DBLINKS
        record: Full gene record (optional, for fallback to unique-id)

    Returns:
        Xref: Primary external reference or None
    """
    # 1. Entrez Gene - Most stable, comprehensive gene database
    # 2. HGNC - Official human gene nomenclature
    # 3. Ensembl - Genomic context
    # 4. UniProt - Protein-centric
    # 5. BioCyc - Pathway database context
    # 6. Other organism-specific databases
    priority_dbs = [
        'ENTREZ', 'NCBI-GENE', 'ENTREZ-GENE',  # Entrez Gene (highest priority)
        'HGNC', 'HGNC-SYMBOL',  # HGNC for human genes
        'ENSEMBL', 'ENSEMBL-PLANT', 'ENSEMBL-PLANTS',  # Ensembl
        'UNIPROT', 'SWISSPROT', 'TREMBL',  # UniProt
        'BIOCYC', 'CHLAMYCYC1',  # BioCyc databases
        'TAIR',  # Arabidopsis
        'MGI', 'SGD', 'FLYBASE', 'ZFIN', 'WORMBASE',  # Model organisms
        'MAIZEGDB', 'GRAMENE-GENES-DB', 'GRAMENE', 'IRGSP', 'PLANTGDB', 'PHYTOZOME', 'TIGR',  # Plants
        'REFSEQ', 'GENBANK', 'EMBL',  # Sequence databases
    ]

    db_mapping = {
        # Gene databases (highest priority)
        'ENTREZ': 'Entrez Gene',
        'NCBI-GENE': 'Entrez Gene',
        'ENTREZ-GENE': 'Entrez Gene',

        # Human gene nomenclature
        'HGNC': 'hgnc',
        'HGNC-SYMBOL': 'hgnc.symbol',

        # Ensembl
        'ENSEMBL': 'ensembl',
        'ENSEMBL-PLANT': 'Ensembl Plants',
        'ENSEMBL-PLANTS': 'Ensembl Plants',

        # UniProt
        'UNIPROT': 'uniprot',
        'SWISSPROT': 'uniprot',
        'TREMBL': 'uniprot',

        # BioCyc
        'BIOCYC': 'biocyc',
        'CHLAMYCYC1': 'biocyc',

        # Plant-specific databases
        'TAIR': 'TAIR gene name',
        'MAIZEGDB': 'Gramene Maize',
        'GRAMENE-GENES-DB': 'Gramene Genes DB',
        'GRAMENE': 'Gramene Genes DB',
        'IRGSP': 'IRGSP Gene',
        'PLANTGDB': 'PlantGDB',
        'PHYTOZOME': 'Phytozome',
        'TIGR': 'TIGR',

        # Model organism databases
        'MGI': 'MGI',
        'SGD': 'SGD',
        'FLYBASE': 'FlyBase',
        'ZFIN': 'ZFIN',
        'WORMBASE': 'WormBase',

        # Sequence databases
        'GENBANK': 'ncbiprotein',
        'REFSEQ': 'refseq',
        'EMBL': 'ena.embl',

        # Other databases
        'KEGG': 'kegg.genes',
        'BRENDA': 'brenda',
        'EC-CODE': 'ec-code',
        'INTERPRO': 'interpro',
        'PFAM': 'pfam',
        'SMART': 'smart',
        'HOMOLOGENE': 'homologene',
        'STRING': 'string'
    }

    # Return the highest priority xref from DBLINKS
    for db_key in priority_dbs:
        if db_key in db_refs:
            return Xref(identifier=db_refs[db_key], dataSource=db_mapping[db_key])

    # If no priority dbs found, return the first available xref from db_mapping
    for db_key, db_id in db_refs.items():
        if db_key in db_mapping:
            return Xref(identifier=db_id, dataSource=db_mapping[db_key])

    # Fallback: Check if UNIQUE-ID matches AT*G pattern for TAIR
    if record and 'UNIQUE-ID' in record:
        unique_id = str(record['UNIQUE-ID'])
        if unique_id.startswith('AT') and 'G' in unique_id:
            return Xref(identifier=unique_id, dataSource='TAIR gene name')

    # Final fallback: use PlantCyc unique-id if available
    if record and 'UNIQUE-ID' in record:
        unique_id = str(record['UNIQUE-ID'])
        return Xref(identifier=unique_id, dataSource='plantcyc')

    return None


def create_gene_properties_from_record(record):
    """
    Create Property elements from gene record data

    Args:
        record: Gene record dictionary

    Returns:
        list: List of Property objects
    """
    properties = []

    if 'UNIQUE-ID' in record:
        properties.append(Property(key='UniqueID', value=str(record['UNIQUE-ID'])))

    # Add ALL database IDs from DBLINKS as individual properties
    dblinks = record.get('DBLINKS', [])
    db_refs = parse_gene_dblinks(dblinks)
    for db_name, db_id in db_refs.items():
        # Create property key from database name
        prop_key = db_name.replace('-', '_').title().replace('_', '')
        properties.append(Property(key=f'AlternativeId_{prop_key}', value=str(db_id)))

    if 'SYNONYMS' in record:
        synonyms = record['SYNONYMS']
        if isinstance(synonyms, list):
            for i, synonym in enumerate(synonyms):
                properties.append(Property(key=f'Synonym_{i + 1}', value=str(synonym).strip()))
        elif isinstance(synonyms, str):
            properties.append(Property(key='Synonym_1', value=synonyms.strip()))

    # Gene-specific properties
    if 'PRODUCT' in record:
        products = record['PRODUCT']
        if isinstance(products, list):
            product_list = ', '.join(str(p) for p in products)
            properties.append(Property(key='Products', value=product_list))
            properties.append(Property(key='ProductCount', value=str(len(products))))
        else:
            properties.append(Property(key='Product', value=str(products)))

    if 'LEFT-END-POSITION' in record:
        properties.append(Property(key='StartPosition', value=str(record['LEFT-END-POSITION'])))

    if 'RIGHT-END-POSITION' in record:
        properties.append(Property(key='EndPosition', value=str(record['RIGHT-END-POSITION'])))

    if 'TRANSCRIPTION-DIRECTION' in record:
        direction = str(record['TRANSCRIPTION-DIRECTION'])
        # Convert to more readable format
        if direction == '+':
            direction = 'Forward'
        elif direction == '-':
            direction = 'Reverse'
        properties.append(Property(key='TranscriptionDirection', value=direction))

    if 'COMPONENT-OF' in record:
        component_of = record['COMPONENT-OF']
        if isinstance(component_of, list):
            component_list = ', '.join(str(c) for c in component_of)
            properties.append(Property(key='ComponentOf', value=component_list))
        else:
            properties.append(Property(key='ComponentOf', value=str(component_of)))

    # Gene location information
    if 'CENTISOME-POSITION' in record:
        properties.append(Property(key='CentisomePosition', value=str(record['CENTISOME-POSITION'])))

    # Gene type information
    if 'TYPES' in record:
        gene_types = record['TYPES']
        if isinstance(gene_types, list):
            types_list = ', '.join(str(t) for t in gene_types)
            properties.append(Property(key='GeneTypes', value=types_list))
        else:
            properties.append(Property(key='GeneType', value=str(gene_types)))

    # Product string
    if 'PRODUCT-STRING' in record:
        properties.append(Property(key='ProductDescription', value=str(record['PRODUCT-STRING'])))

    # Operon inf
    if 'COMPONENT-OF' in record:
        operon_info = record['COMPONENT-OF']
        if isinstance(operon_info, list):
            operon_list = ', '.join(str(o) for o in operon_info)
            properties.append(Property(key='Operons', value=operon_list))
        else:
            properties.append(Property(key='Operon', value=str(operon_info)))

    # Additional gene feats
    if 'INTERRUPTED?' in record:
        properties.append(Property(key='Interrupted', value=str(record['INTERRUPTED?'])))

    # ========================================================================
    # Dynamic parsing for ALL other fields not explicitly handled above
    # ========================================================================
    skip_fields = {
        'UNIQUE-ID', 'DBLINKS', 'SYNONYMS', 'PRODUCT', 'LEFT-END-POSITION',
        'RIGHT-END-POSITION', 'TRANSCRIPTION-DIRECTION', 'COMPONENT-OF',
        'CENTISOME-POSITION', 'TYPES', 'PRODUCT-STRING', 'INTERRUPTED?',
        'COMMON-NAME', 'CITATIONS', 'COMMENT', 'CREDITS',
    }
    dynamic_props = create_dynamic_properties(record, skip_fields=skip_fields)
    properties.extend(dynamic_props)

    return properties


def create_gene_comments_from_record(record):
    """
    Create Comment elements from gene record data

    Args:
        record: Gene record dictionary

    Returns:
        list: List of Comment objects
    """
    comments = []

    #get source
    sources = record.get('CREDITS', [])
    source_str = ", ".join(str(s) for s in sources) if isinstance(sources, list) else str(sources)

    # Add main COMMENT field (descriptive text only)
    if 'COMMENT' in record:
        comment_text = record['COMMENT']
        if isinstance(comment_text, list):
            for c in comment_text:
                # Clean CITS tags from comment text
                cleaned_text = re.sub(r'\|CITS:\s*\[(\d+)\]\|', '', str(c)).strip()
                if cleaned_text:
                    comments.append(Comment(value=cleaned_text, source=source_str))
        elif isinstance(comment_text, str):
            # Clean CITS tags from comment text
            cleaned_text = re.sub(r'\|CITS:\s*\[(\d+)\]\|', '', comment_text).strip()
            if cleaned_text:
                comments.append(Comment(value=cleaned_text, source=source_str))

    return comments


def create_gene_citations_from_record(record):
    """
    Create Citation and CitationRef elements from gene record

    Args:
        record: Gene record dictionary

    Returns:
        tuple: (citations_list, citation_refs_list)
    """
    citations = []
    citation_refs = []

    if 'CITATIONS' in record:
        citation_data = record['CITATIONS']
        if isinstance(citation_data, str):
            citation_data = [citation_data]
        elif not isinstance(citation_data, list):
            citation_data = []

        for citation_entry in citation_data:
            citation_parts = str(citation_entry).split(':')
            if len(citation_parts) >= 1:
                citation_id = citation_parts[0]

                citation = Citation(
                    elementId=f"citation_{citation_id}",
                    xref=Xref(identifier=citation_id, dataSource="PubMed")
                )
                citations.append(citation)
                citation_ref = CitationRef(elementRef=f"citation_{citation_id}")
                citation_refs.append(citation_ref)

    return citations, citation_refs


def create_enhanced_datanode_from_gene(record, citation_manager=None):
    """
    Convert a gene record to a comprehensive GPML DataNode

    Args:
        record: Dictionary containing gene data from parsing_utils
        citation_manager: CitationManager instance (optional)

    Returns:
        tuple: (DataNode, citations list)
    """
    element_id = record.get('UNIQUE-ID', str(uuid.uuid4()))
    raw_label = record.get('COMMON-NAME', record.get('UNIQUE-ID', 'Unknown Compound'))
    text_label = clean_text_label(raw_label)

    # Parse db links
    dblinks = record.get('DBLINKS', [])
    db_refs = parse_gene_dblinks(dblinks)
    primary_xref = get_primary_gene_xref(db_refs, record)

    # Create props and citations
    properties = create_gene_properties_from_record(record)
    comments = create_gene_comments_from_record(record)

    # Handle citations
    if citation_manager:
        citation_refs = create_citation_refs_from_record(record, citation_manager)
    else:
        citations, citation_refs = create_gene_citations_from_record(record)

    # Handle annotations (Species)
    annotations, annotation_refs = create_species_annotation(record)

    # Use  gene graphics
    graphics = standard_graphics.create_gene_graphics(0.0, 0.0)

    datanode = DataNode(
        elementId=element_id,
        textLabel=str(text_label),
        type=DataNodeType.GENE_PRODUCT,
        xref=primary_xref,
        states=[],
        graphics=graphics,
        comments=comments,
        properties=properties,
        citationRefs=citation_refs if citation_manager else citation_refs,
        annotationRefs=annotation_refs
    )

    return (datanode, [], annotations) if citation_manager else (datanode, citations, annotations)


def create_citation_refs_from_record(record, citation_manager):
    """
    Create CitationRef elements from gene record using CitationManager.
    Extracts citations from both CITATIONS field and inline citations in COMMENT field.

    Args:
        record: Gene record dictionary
        citation_manager: CitationManager instance

    Returns:
        list: List of CitationRef objects (deduplicated)
    """
    if not citation_manager:
        return []

    all_citations = []
    gene_id = record.get('UNIQUE-ID', '')

    # 1. Extract citations from CITATIONS field
    if 'CITATIONS' in record:
        citations = record['CITATIONS']
        if not isinstance(citations, list):
            citations = [citations] if citations else []
        all_citations.extend([str(c) for c in citations if c])

    # 2. Extract inline citations from COMMENT field
    # Pattern: |CITS:[12345]| or |CITS: [12345]|
    if 'COMMENT' in record:
        comment_text = record['COMMENT']
        if isinstance(comment_text, list):
            comment_text = ' '.join(str(c) for c in comment_text)
        else:
            comment_text = str(comment_text)

        # Find all inline citations
        inline_pattern = r'\|CITS:\s*\[(\d+)\]\|'
        inline_citations = re.findall(inline_pattern, comment_text)
        all_citations.extend(inline_citations)

    # 3. Deduplicate and process through CitationManager
    unique_citations = list(dict.fromkeys(all_citations))  # Preserves order while deduplicating

    if unique_citations:
        citation_refs = citation_manager.process_element_citations(gene_id, unique_citations)
        return citation_refs

    return []


def create_enhanced_datanodes_from_genes(genes_file, citation_manager=None):
    """
    Create comprehensive DataNodes from a genes file with citation support.

    Args:
        genes_file: Path to the genes.dat file
        citation_manager: CitationManager instance (optional)

    Returns:
        tuple: (list of DataNode objects, list of all Citation objects)
    """
    processor = parsing_utils.read_and_parse(genes_file)

    datanodes = []
    all_citations = []
    all_annotations = []

    # Statistics tracking - count per unique record
    xref_stats = {
        'total_unique_records': 0,
        'records_with_xref': 0,
        'xref_by_database': {}
    }

    citation_stats = {
        'records_with_citations': 0,
        'total_citation_refs': 0,
        'unique_pubmed_ids': set()
    }

    for record in processor.records:
        xref_stats['total_unique_records'] += 1

        # Get primary xref for statistics
        dblinks = record.get('DBLINKS', [])
        db_refs = parse_gene_dblinks(dblinks)
        primary_xref = get_primary_gene_xref(db_refs, record)

        if primary_xref:
            xref_stats['records_with_xref'] += 1
            db_source = primary_xref.dataSource
            xref_stats['xref_by_database'][db_source] = xref_stats['xref_by_database'].get(db_source, 0) + 1

        # Track citations from both CITATIONS field and inline citations
        all_citations = []

        # 1. Get citations from CITATIONS field
        if 'CITATIONS' in record:
            citations = record['CITATIONS']
            if not isinstance(citations, list):
                citations = [citations] if citations else []
            all_citations.extend([str(c) for c in citations if c])

        # 2. Extract inline citations from COMMENT field
        if 'COMMENT' in record:
            comment_text = record['COMMENT']
            if isinstance(comment_text, list):
                comment_text = ' '.join(str(c) for c in comment_text)
            else:
                comment_text = str(comment_text)

            # Find all inline citations using same pattern as create_citation_refs_from_record
            inline_pattern = r'\|CITS:\s*\[(\d+)\]\|'
            inline_citations = re.findall(inline_pattern, comment_text)
            all_citations.extend(inline_citations)

        # Deduplicate and track stats
        unique_citations = list(dict.fromkeys(all_citations))

        if unique_citations:
            citation_stats['records_with_citations'] += 1
            citation_stats['total_citation_refs'] += len(unique_citations)

            # Extract PubMed IDs
            for citation in unique_citations:
                citation_str = str(citation)
                if citation_str.isdigit():
                    citation_stats['unique_pubmed_ids'].add(citation_str)

        # Create datanode
        datanode, citations, annotations = create_enhanced_datanode_from_gene(record, citation_manager)
        datanodes.append(datanode)
        if citations:
            all_citations.extend(citations)
        if annotations:
            all_annotations.extend(annotations)

    # Print statistics
    print("\n" + "="*60)
    print("GENE XREF ANNOTATION STATISTICS")
    print("="*60)
    print(f"Total unique gene records: {xref_stats['total_unique_records']}")
    print(f"Records with xrefs: {xref_stats['records_with_xref']}")
    if xref_stats['total_unique_records'] > 0:
        print(f"Annotation rate: {xref_stats['records_with_xref']/xref_stats['total_unique_records']*100:.1f}%")
    print(f"\nPrimary xrefs by database (all types):")
    for db_name in sorted(xref_stats['xref_by_database'].keys()):
        count = xref_stats['xref_by_database'][db_name]
        percentage = (count / xref_stats['total_unique_records'] * 100) if xref_stats['total_unique_records'] > 0 else 0
        print(f"  {db_name}: {count} ({percentage:.1f}%)")

    print(f"\nCitation Statistics:")
    print(f"Records with citations: {citation_stats['records_with_citations']} ({citation_stats['records_with_citations']/xref_stats['total_unique_records']*100:.1f}%)")
    print(f"Total citation references: {citation_stats['total_citation_refs']}")
    print(f"Unique PubMed IDs: {len(citation_stats['unique_pubmed_ids'])}")
    print("="*60 + "\n")

    return datanodes, all_citations, all_annotations


def print_gene_datanode_summary(datanode: DataNode, index):
    """
    Print a summary of a gene DataNode

    Args:
        datanode: DataNode object
        index: Index number for display
    """
    print(f"\nGene DataNode {index + 1}:")
    print(f"  ID: {datanode.elementId}")
    print(f"  Label: {datanode.textLabel}")
    print(f"  Type: {datanode.type.value}")

    if datanode.xref:
        print(f"  Primary Xref: {datanode.xref.dataSource}:{datanode.xref.identifier}")
    else:
        print(f"  Primary Xref: None")

    print(f"  Properties: {len(datanode.properties)} items")
    for prop in datanode.properties[:3]:
        print(f"    - {prop.key}: {prop.value}")
    if len(datanode.properties) > 3:
        print(f"    ... and {len(datanode.properties) - 3} more")

    print(f"  Comments: {len(datanode.comments)} items")
    if datanode.comments:
        for comment in datanode.comments:
            comment_text = comment.value[:50] + '...' if len(comment.value) > 50 else comment.value
            print(f"    - {comment.source}: {comment_text}")

    if datanode.citationRefs:
        print(f"  Citations: {len(datanode.citationRefs)} references")


if __name__ == "__main__":
    datanodes, citations, annotations = create_enhanced_datanodes_from_genes("genes.dat")
    print(f"Created {len(datanodes)} enhanced Gene DataNodes")
    print(f"Found {len(citations)} citations")
    print(f"Found {len(annotations)} annotations")

    for i, node in enumerate(datanodes[:3]):
        print_gene_datanode_summary(node, i)

    print("\n" + "=" * 80)
    if datanodes:
        gpml_writer = GPMLWriter()
        sample_node = datanodes[0]
        gpml_output = gpml_writer.write_datanode(sample_node)
        print("Sample gene GPML output:")
        print(gpml_output[:500] if len(gpml_output) > 500 else gpml_output)