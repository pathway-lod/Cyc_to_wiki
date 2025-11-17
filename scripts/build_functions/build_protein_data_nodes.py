from scripts.data_structure.wiki_data_structure import (
    Xref, Graphics, DataNode, DataNodeType, HAlign, VAlign,
    BorderStyle, ShapeType, Property, Comment, CitationRef, Citation, Group, GroupType
)
from scripts.parsing_functions import parsing_utils
from scripts.utils.HTML_cleaner import clean_text_label
from scripts.utils import standard_graphics
from scripts.utils.property_parser import (
    create_dynamic_properties, handle_unique_id, handle_synonyms,
    handle_go_terms, create_numbered_properties
)
from scripts.object2gmpl.gpml_writer import GPMLWriter
import uuid
import re


def parse_protein_dblinks(dblinks):
    """
    Parse DBLINKS to extract database references for proteins

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


def get_primary_protein_xref(db_refs, record=None):
    """
    Get the primary external reference for a protein based on BridgeDb protein databases

    Args:
        db_refs: Dictionary of database references
        record: Full protein record (optional, for fallback to unique-id)

    Returns:
        Xref: Primary external reference or None
    """
    # BridgeDb protein database mapping
    # Priority order based on WikiPathways best practices:
    # 1. Entrez Gene - Most stable, comprehensive gene database
    # 2. HGNC - Official human gene nomenclature
    # 3. Ensembl - Genomic context
    # 4. UniProt - Protein-centric
    # 5. BioCyc - Pathway database context
    # 6. Other databases
    priority_dbs = [
        'ENTREZ', 'NCBI-GENE', 'ENTREZ-GENE',  # Entrez Gene (highest priority)
        'HGNC', 'HGNC-SYMBOL',  # HGNC for human genes (for later versions of other cycs)
        'ENSEMBL', 'ENSEMBL-PLANT', 'ENSEMBL-PLANTS',  # Ensembl
        'UNIPROT', 'SWISS-PROT', 'TREMBL',  # UniProt variants
        'BIOCYC', 'CHLAMYCYC1',  # BioCyc databases
        'NCBI-PROTEIN-ID', 'NCBI-PROTEIN', 'REFSEQ',  # Other NCBI databases
        'PDB', 'INTERPRO', 'PFAM',  # Structure and domain databases
        'IPI', 'BRENDA', 'MACIE', 'RESID', 'SMART', 'SUPFAM',  # Other protein databases
        'PHOSPHOSITE', 'STRING', 'SUBSTRATEB', 'TTD',
        'SWISS-MODEL', 'SPRINT',
        'GENBANK', 'PIR',
        'TAIR', 'TAIR-PROTEIN', 'ARAPORT',  # Plant-specific databases (ARAPORT -> TAIR)
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

        # UniProt variants
        'UNIPROT': 'uniprot',
        'SWISS-PROT': 'uniprot',
        'TREMBL': 'uniprot',

        # BioCyc
        'BIOCYC': 'biocyc',
        'CHLAMYCYC1': 'biocyc',

        # BridgeDb protein resources
        'BRENDA': 'brenda',
        'INTERPRO': 'interpro',
        'IPI': 'IPI',
        'MACIE': 'macie',
        'NCBI-PROTEIN-ID': 'ncbiprotein',
        'NCBI-PROTEIN': 'ncbiprotein',
        'PDB': 'pdb',
        'PFAM': 'pfam',
        'PHOSPHOSITE': 'phosphosite.protein',
        'RESID': 'resid',
        'SMART': 'smart',
        'SPRINT': 'SPRINT',
        'STRING': 'string',
        'SUBSTRATEB': 'SubstrateDB',
        'SUPFAM': 'supfam',
        'SWISS-MODEL': 'swiss-model',
        'TTD': 'ttd.target',

        # Plant-specific protein databases
        'TAIR': 'TAIR',
        'TAIR-PROTEIN': 'TAIR',
        'ARAPORT': 'TAIR',

        # Other sequence databases
        'REFSEQ': 'refseq',
        'GENBANK': 'GenBank',
        'PIR': 'PIR'
    }

    # Return the highest priority xref
    for db_key in priority_dbs:
        if db_key in db_refs:
            if db_key in db_mapping:
                return Xref(identifier=db_refs[db_key], dataSource=db_mapping[db_key])

    # If no priority dbs found, return the first available xref
    for db_key, db_id in db_refs.items():
        if db_key in db_mapping:
            return Xref(identifier=db_id, dataSource=db_mapping[db_key])

    # Final fallback: use PlantCyc unique-id if available
    if record and 'UNIQUE-ID' in record:
        unique_id = str(record['UNIQUE-ID'])
        return Xref(identifier=unique_id, dataSource='plantcyc')

    return None


def is_protein_complex(record):
    """
    Determine if a protein record represents a complex.

    Args:
        record: Protein record dictionary

    Returns:
        bool: True if this is a protein complex
    """
    # First check: Does it have COMPONENTS field? This is the most reliable indicator
    if 'COMPONENTS' in record:
        components = record.get('COMPONENTS', [])
        if components:  # Has non-empty components
            return True

    # Second check: Check TYPES field for complex-related keywords
    types = record.get('TYPES', [])
    if not isinstance(types, list):
        types = [types] if types else []

    # Check if any type indicates this is a complex
    for type_val in types:
        if isinstance(type_val, str):
            type_str = type_val.upper()
            if 'PROTEIN-COMPLEXES' in type_str:
                return True

    return False


def parse_complex_components(record):
    """
    Parse the COMPONENTS field to get monomer IDs.

    Args:
        record: Protein complex record dictionary

    Returns:
        list: List of component monomer IDs
    """
    components = record.get('COMPONENTS', [])
    if not isinstance(components, list):
        components = [components] if components else []

    # Clean component IDs
    component_ids = []
    for comp in components:
        if comp:
            component_ids.append(str(comp).strip())

    return component_ids


def create_protein_properties_from_record(record):
    """
    Create Property elements from protein record data:
    - Custom handlers for special fields (DBLINKS, SYNONYMS, COMPONENTS, CATALYZES)
    - Dynamic parsing for all other fields
    """
    properties = []

    # ========================================================================
    # PART 1: Critical fields with custom handling
    # ========================================================================

    # 1. UNIQUE-ID
    if 'UNIQUE-ID' in record:
        properties.extend(handle_unique_id(record['UNIQUE-ID'], record))

    # 2. DBLINKS - Handle specially for EC numbers and alternative IDs
    dblinks = record.get('DBLINKS', [])
    db_refs = parse_protein_dblinks(dblinks)
    if 'EC-CODE' in db_refs or 'BRENDA' in db_refs:
        ec_number = db_refs.get('EC-CODE') or db_refs.get('BRENDA')
        properties.append(Property(key='ECNumber', value=str(ec_number)))

    for db_name, db_id in db_refs.items():
        if db_name not in ['EC-CODE', 'BRENDA']:
            prop_key = db_name.replace('-', '_').title().replace('_', '')
            properties.append(Property(key=f'AlternativeId_{prop_key}', value=str(db_id)))

    # 3. SYNONYMS - Create numbered properties
    if 'SYNONYMS' in record:
        properties.extend(handle_synonyms(record['SYNONYMS'], record))

    # 4. COMPONENTS - Special handling for complexes
    if is_protein_complex(record):
        if 'COMPONENTS' in record:
            components = parse_complex_components(record)
            if components:
                properties.append(Property(key='ComplexComponents', value=', '.join(components)))
                properties.append(Property(key='ComponentCount', value=str(len(components))))
        properties.append(Property(key='IsComplex', value='true'))

    # 5. CATALYZES - Special handling to clean ENZ prefixes
    if 'CATALYZES' in record:
        catalyzes = record['CATALYZES']

        def remove_enz_prefix(reaction_id):
            """Remove ENZ prefix from reaction IDs to match pathway reactions"""
            reaction_id = str(reaction_id).strip()
            if reaction_id.startswith('ENZRXN-'):
                return reaction_id[7:]
            elif reaction_id.startswith('ENZRXNQT-'):
                return reaction_id[9:]
            elif reaction_id.startswith('ENZRXNIO2-'):
                return reaction_id[10:]
            elif reaction_id.startswith('ENZRXN'):
                return reaction_id[6:]
            else:
                return reaction_id

        if isinstance(catalyzes, list):
            cleaned_reactions = []
            for rxn_id in catalyzes:
                cleaned_id = remove_enz_prefix(rxn_id)
                if not cleaned_id.startswith('RXN'):
                    cleaned_id = f'RXN-{cleaned_id}'
                cleaned_reactions.append(cleaned_id)

            catalyzes_list = ', '.join(cleaned_reactions)
            properties.append(Property(key='Catalyzes', value=catalyzes_list))
            properties.append(Property(key='CatalyzesCount', value=str(len(catalyzes))))
        else:
            cleaned_id = remove_enz_prefix(catalyzes)
            if not cleaned_id.startswith('RXN'):
                cleaned_id = f'RXN-{cleaned_id}'
            properties.append(Property(key='Catalyzes', value=cleaned_id))

    # 6. GO-TERMS - Keep explicit for clarity
    if 'GO-TERMS' in record:
        properties.extend(handle_go_terms(record['GO-TERMS'], record))

    # ========================================================================
    # PART 2: Dynamic parsing for ALL other fields
    # ========================================================================

    # Fields to skip
    skip_fields = {
        'UNIQUE-ID', 'DBLINKS', 'SYNONYMS', 'COMPONENTS', 'CATALYZES', 'GO-TERMS',
        'COMMON-NAME',  # Already in textLabel
        'CITATIONS',    # Handled separately
        'COMMENT',      # Handled separately
        'CREDITS',      # Handled in comment source
    }

    # Parse all remaining fields dynamically
    dynamic_props = create_dynamic_properties(record, skip_fields=skip_fields)
    properties.extend(dynamic_props)

    return properties


def create_complex_group(record, citation_manager=None):
    """
    Create a Group object for a protein complex.

    Args:
        record: Protein complex record dictionary
        citation_manager: CitationManager instance (optional)

    Returns:
        Group: Group object representing the complex
    """
    element_id = record.get('UNIQUE-ID', str(uuid.uuid4()))
    raw_label = record.get('COMMON-NAME', record.get('UNIQUE-ID', 'Protein Complex'))
    text_label = clean_text_label(raw_label)

    # Parse database links
    dblinks = record.get('DBLINKS', [])
    db_refs = parse_protein_dblinks(dblinks)
    primary_xref = get_primary_protein_xref(db_refs, record)

    # Create properties and comments
    properties = create_protein_properties_from_record(record)
    comments = create_protein_comments_from_record(record)

    # Handle citations
    citation_refs = []
    if citation_manager:
        citation_refs = create_citation_refs_from_record(record, citation_manager)

    # Create graphics for the group
    graphics = standard_graphics.create_complex_group_graphics(0.0, 0.0, 200.0, 100.0)

    group = Group(
        elementId=element_id,
        textLabel=str(text_label),
        type=GroupType.COMPLEX,
        xref=primary_xref,
        graphics=graphics,
        comments=comments,
        properties=properties,
        citationRefs=citation_refs
    )

    return group


def create_monomer_datanode(monomer_id, parent_complex_id, monomer_records=None):
    """
    Create a DataNode for a monomer that's part of a complex.

    Args:
        monomer_id: ID of the monomer
        parent_complex_id: ID of the parent complex (for groupRef)
        monomer_records: Optional dict of all protein records for lookup

    Returns:
        DataNode: DataNode for the monomer
    """
    # Try to find the monomer record if available
    monomer_record = None
    if monomer_records and monomer_id in monomer_records:
        monomer_record = monomer_records[monomer_id]

    # Create basic datanode
    element_id = monomer_id
    text_label = monomer_id  # Default to ID

    if monomer_record:
        # Use actual record data if available
        raw_label = monomer_record.get('COMMON-NAME', monomer_record.get('UNIQUE-ID', monomer_id))
        text_label = clean_text_label(raw_label)

        # Get xref
        dblinks = monomer_record.get('DBLINKS', [])
        db_refs = parse_protein_dblinks(dblinks)
        primary_xref = get_primary_protein_xref(db_refs, monomer_record)

        properties = create_protein_properties_from_record(monomer_record)
        comments = create_protein_comments_from_record(monomer_record)
    else:
        # Minimal properties for unknown monomers
        primary_xref = None
        properties = [Property(key='UniqueID', value=monomer_id)]
        comments = []

    # Use protein graphics
    graphics = standard_graphics.create_protein_graphics(0.0, 0.0)

    datanode = DataNode(
        elementId=element_id,
        textLabel=str(text_label),
        type=DataNodeType.PROTEIN,
        groupRef=parent_complex_id,  # Reference to parent complex
        xref=primary_xref,
        states=[],
        graphics=graphics,
        comments=comments,
        properties=properties,
        citationRefs=[]
    )

    return datanode
def create_protein_comments_from_record(record):
    """
    Create Comment elements from protein record data

    Args:
        record: Protein record dictionary

    Returns:
        list: List of Comment objects
    """
    comments = []

    # Try to get source information
    sources = record.get('CREDITS', [])
    source_str = ", ".join(str(s) for s in sources) if isinstance(sources, list) else str(sources)

    # Add main COMMENT field (descriptive text only)
    if 'COMMENT' in record:
        comment_text = record['COMMENT']
        if isinstance(comment_text, list):
            for c in comment_text:
                comments.append(Comment(value=str(c).strip(), source=source_str))
        elif isinstance(comment_text, str):
            comments.append(Comment(value=comment_text.strip(), source=source_str))

    return comments


def create_protein_citations_from_record(record):
    """
    Create Citation and CitationRef elements from protein record

    Args:
        record: Protein record dictionary

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


def create_citation_refs_from_record(record, citation_manager):
    """
    Create CitationRef elements from protein record using CitationManager.
    Extracts citations from both CITATIONS field and inline citations in COMMENT field.

    Args:
        record: Protein record dictionary
        citation_manager: CitationManager instance

    Returns:
        list: List of CitationRef objects (deduplicated)
    """
    if not citation_manager:
        return []

    all_citations = []
    protein_id = record.get('UNIQUE-ID', '')

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
        citation_refs = citation_manager.process_element_citations(protein_id, unique_citations)
        return citation_refs

    return []


def create_enhanced_datanode_from_protein(record, citation_manager=None):
    """Original function for creating protein DataNodes """
    element_id = record.get('UNIQUE-ID', str(uuid.uuid4()))
    raw_label = record.get('COMMON-NAME', record.get('UNIQUE-ID', 'Unknown Protein'))
    text_label = clean_text_label(raw_label)

    dblinks = record.get('DBLINKS', [])
    db_refs = parse_protein_dblinks(dblinks)
    primary_xref = get_primary_protein_xref(db_refs, record)

    properties = create_protein_properties_from_record(record)
    comments = create_protein_comments_from_record(record)

    citation_refs = []
    if citation_manager:
        citation_refs = create_citation_refs_from_record(record, citation_manager)

    # Use standard protein graphics
    graphics = standard_graphics.create_protein_graphics(0.0, 0.0)

    datanode = DataNode(
        elementId=element_id,
        textLabel=str(text_label),
        type=DataNodeType.PROTEIN,
        xref=primary_xref,
        states=[],
        graphics=graphics,
        comments=comments,
        properties=properties,
        citationRefs=citation_refs
    )

    return (datanode, [])

def create_enhanced_datanodes_from_proteins(proteins_file, citation_manager=None):
    """
    Create comprehensive DataNodes and Groups from a proteins file with complex support.

    Args:
        proteins_file: Path to the proteins.dat file
        citation_manager: CitationManager instance (optional)

    Returns:
        tuple: (list of DataNode objects, list of Group objects, list of all Citation objects)
    """
    processor = parsing_utils.read_and_parse(proteins_file)

    datanodes = []
    groups = []
    all_citations = []

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

    # Track unique species/taxon
    unique_species = set()

    # First pass: build index of all protein records and collect stats
    protein_records = {}
    for record in processor.records:
        unique_id = record.get('UNIQUE-ID')
        if unique_id:
            protein_records[unique_id] = record
            xref_stats['total_unique_records'] += 1

            # Get primary xref for this unique record
            dblinks = record.get('DBLINKS', [])
            db_refs = parse_protein_dblinks(dblinks)
            primary_xref = get_primary_protein_xref(db_refs, record)

            if primary_xref:
                xref_stats['records_with_xref'] += 1
                db_source = primary_xref.dataSource
                xref_stats['xref_by_database'][db_source] = xref_stats['xref_by_database'].get(db_source, 0) + 1

            # Track species
            if 'SPECIES' in record:
                species = record['SPECIES']
                if isinstance(species, list):
                    for s in species:
                        unique_species.add(str(s))
                else:
                    unique_species.add(str(species))

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
                    # PubMed IDs are typically just numbers
                    if citation_str.isdigit():
                        citation_stats['unique_pubmed_ids'].add(citation_str)

    # Second pass: process proteins and complexes (create actual DataNodes/Groups)
    processed_monomers = set()  # Track monomer IDs we've already added as part of complexes
    processed_protein_ids = set()  # Track protein unique IDs we've already created nodes for


    for unique_id, record in protein_records.items():
        if unique_id in processed_protein_ids:
            continue
        processed_protein_ids.add(unique_id)

        if is_protein_complex(record):
            # Create a Group for the complex
            complex_group = create_complex_group(record, citation_manager)
            groups.append(complex_group)

            # Create DataNodes for each component
            component_ids = parse_complex_components(record)
            for comp_id in component_ids:
                if comp_id not in processed_monomers:
                    monomer_node = create_monomer_datanode(
                        comp_id,
                        complex_group.elementId,
                        protein_records
                    )
                    datanodes.append(monomer_node)
                    processed_monomers.add(comp_id)
        else:
            # Regular protein - create as DataNode (only if not already part of a complex)
            if unique_id not in processed_monomers:
                datanode = create_enhanced_datanode_from_protein(record, citation_manager)
                protein_node = datanode[0] if isinstance(datanode, tuple) else datanode
                datanodes.append(protein_node)

    # Final deduplication pass: keep only first occurrence of each elementId
    # (some monomers may be components of multiple complexes)
    seen_element_ids = set()
    deduplicated_datanodes = []
    for datanode in datanodes:
        if datanode.elementId not in seen_element_ids:
            deduplicated_datanodes.append(datanode)
            seen_element_ids.add(datanode.elementId)

    if len(deduplicated_datanodes) < len(datanodes):
        print(f"\nDeduplication: Removed {len(datanodes) - len(deduplicated_datanodes)} duplicate DataNodes")
        datanodes = deduplicated_datanodes

    # Print statistics
    print("\n" + "="*60)
    print("PROTEIN XREF ANNOTATION STATISTICS")
    print("="*60)
    print(f"Total unique protein records: {xref_stats['total_unique_records']}")
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

    print(f"\nSpecies/Taxon Information:")
    print(f"Unique species found: {len(unique_species)}")
    print("="*60 + "\n")

    return datanodes, groups, all_citations




def print_protein_datanode_summary(datanode: DataNode, index):
    """
    Print a summary of a protein DataNode

    Args:
        datanode: DataNode object
        index: Index number for display
    """
    print(f"\nProtein DataNode {index + 1}:")
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
    datanodes, citations = create_enhanced_datanodes_from_proteins("proteins.dat")
    print(f"Created {len(datanodes)} enhanced Protein DataNodes")
    print(f"Found {len(citations)} citations")

    for i, node in enumerate(datanodes[:3]):
        print_protein_datanode_summary(node, i)

    # Test GPML output
    print("\n" + "=" * 80)
    if datanodes:
        gpml_writer = GPMLWriter()
        sample_node = datanodes[0]
        gpml_output = gpml_writer.write_datanode(sample_node)
        print("Sample protein GPML output:")
        print(gpml_output[:] if len(gpml_output) > 500 else gpml_output)