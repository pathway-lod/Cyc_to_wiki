from scripts.data_structure.wiki_data_structure import (
    Xref, Graphics, DataNode, DataNodeType, HAlign, VAlign,
    BorderStyle, ShapeType, Property, Comment, CitationRef
)
from scripts.parsing_functions import parsing_utils
from scripts.utils.HTML_cleaner import clean_text_label
from scripts.utils import standard_graphics
from scripts.object2gmpl.gpml_writer import GPMLWriter
import uuid
import re


def parse_dblinks(dblinks):
    """
    Parse DBLINKS to extract database references

    Args:
        dblinks: String, list of strings, or list of tuples

    Returns:
        dict: Dictionary with database names as keys and IDs as values
    """
    db_refs = {}

    # If dblinks is a string, convert to list
    if isinstance(dblinks, str):
        dblinks = [dblinks]

    # If dblinks is not a list, return empty dict
    if not isinstance(dblinks, list):
        return db_refs

    for dblink in dblinks:
        # Handle both string format and tuple/list format
        if isinstance(dblink, str):
            try:
                # Match database names with hyphens and pipes (e.g., LIGAND-CPD, |Wikipedia|)
                match = re.match(r'\((\|?[A-Za-z0-9-]+\|?)\s+"([^"]+)"', dblink)
                if match:
                    db_name = match.group(1).strip('|').upper()  # Remove pipes and uppercase
                    db_id = match.group(2)
                    db_refs[db_name] = db_id
            except Exception as e:
                print(f"Warning: Could not parse dblink string: {dblink}")
        elif isinstance(dblink, (list, tuple)) and len(dblink) >= 2:
            # Handle parsed format like ('PUBCHEM', '25245269', ...)
            try:
                db_name = str(dblink[0]).upper()
                db_id = str(dblink[1])
                db_refs[db_name] = db_id
            except Exception as e:
                print(f"Warning: Could not parse dblink tuple: {dblink}")
        else:
            # Unknown format
            print(f"DEBUG: Unknown dblink format - type: {type(dblink)}, value: {repr(dblink)[:100]}")

    return db_refs


def get_primary_xref(db_refs, inchi_key=None):
    """
    Get the primary external reference for a compound

    Args:
        db_refs: Dictionary of database references
        inchi_key: InChI key from the record (optional, takes priority)

    Returns:
        Xref: Primary external reference or None
    """
    # If InChI key is available, clean it and use it as primary xref
    if inchi_key:
        # Remove "InChIKey=" prefix if present and clean whitespace
        clean_inchi_key = str(inchi_key).replace('InChIKey=', '').strip()
        # Only return if we have a non-empty key after cleaning
        if clean_inchi_key:
            return Xref(identifier=clean_inchi_key, dataSource='InChIKey')

    priority_dbs = ['PUBCHEM', 'CHEBI', 'HMDB', 'CHEMSPIDER', 'LIGAND-CPD', 'CAS', 'CHEMBL', 'DRUGBANK',
                    'KEGG-DRUG', 'KEGG-GLYCAN', 'LIPIDMAPS', 'LIPIDBANK', 'PHARMGKB',
                    'SWISSLIPIDS', 'TTD', 'WIKIDATA',
                    'METANETX', 'SEED', 'BIGG', 'ZINC',
                    'CHEBI']

    # Mapping from flat file database names to full database names
    db_mapping = {
        'PUBCHEM': 'PubChem-compound',
        'CHEBI': 'ChEBI',
        'HMDB': 'HMDB',
        'CHEMSPIDER': 'ChemSpider',
        'LIGAND-CPD': 'KEGG Compound',
        'CAS': 'CAS',
        'CHEMBL': 'ChEMBL',
        'DRUGBANK': 'DrugBank',
        'KEGG-DRUG': 'KEGG Drug',
        'KEGG-GLYCAN': 'KEGG Glycan',
        'LIPIDMAPS': 'LIPID MAPS',
        'LIPIDBANK': 'LipidBank',
        'PHARMGKB': 'PharmGKB',
        'SWISSLIPIDS': 'SwissLipids',
        'TTD': 'TTD',
        'WIKIDATA': 'Wikidata',
        # Non-BridgeDb databases (fallback)
        'METANETX': 'MetaNetX',
        'SEED': 'ModelSEED',
        'BIGG': 'BiGG',
        'ZINC': 'ZINC'
    }

    # Return the highest priority xref
    for db_key in priority_dbs:
        if db_key in db_refs:
            return Xref(identifier=db_refs[db_key], dataSource=db_mapping[db_key])

    # If no priority dbs found, return the first available xref that's in BridgeDb
    for db_key, db_id in db_refs.items():
        if db_key in db_mapping:
            return Xref(identifier=db_id, dataSource=db_mapping[db_key])

    return None


def create_properties_from_record(record):
    """
    Create Property elements from compound record data

    Args:
        record: Compound record dictionary

    Returns:
        list: List of Property objects
    """
    properties = []

    if 'UNIQUE-ID' in record:
        properties.append(Property(key='UniqueID', value=str(record['UNIQUE-ID'])))

    # Add ALL database IDs from DBLINKS as individual properties
    dblinks = record.get('DBLINKS', [])
    db_refs = parse_dblinks(dblinks)
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

    if 'MOLECULAR-WEIGHT' in record:
        properties.append(Property(key='MolecularWeight', value=str(record['MOLECULAR-WEIGHT'])))

    if 'MONOISOTOPIC-MW' in record:
        properties.append(Property(key='MonoisotopicWeight', value=str(record['MONOISOTOPIC-MW'])))

    if 'CHEMICAL-FORMULA' in record:
        formula_parts = record['CHEMICAL-FORMULA']
        if isinstance(formula_parts, list):
            formula = ''
            for part in formula_parts:
                if isinstance(part, tuple) and len(part) == 2:
                    element, count = part
                    formula += f" {element} {count} "
                else:
                    formula += str(part)
            properties.append(Property(key='ChemicalFormula', value=formula.strip()))

    if 'INCHI' in record:
        properties.append(Property(key='InChI', value=str(record['INCHI'])))

    if 'INCHI-KEY' in record:
        properties.append(Property(key='InChIKey', value=str(record['INCHI-KEY'])))

    if 'SMILES' in record:
        properties.append(Property(key='SMILES', value=str(record['SMILES'])))

    if 'GIBBS-0' in record:
        properties.append(Property(key='GibbsFreeEnergy', value=str(record['GIBBS-0'])))

    if 'TYPES' in record:
        compound_type = record['TYPES']
        if isinstance(compound_type, list):
            compound_type = ', '.join(str(t) for t in compound_type)
        properties.append(Property(key='CompoundType', value=str(compound_type)))

    if 'ABBREV-NAME' in record:
        properties.append(Property(key='AbbreviatedName', value=str(record['ABBREV-NAME'])))

    if 'SYSTEMATIC-NAME' in record:
        properties.append(Property(key='SystematicName', value=str(record['SYSTEMATIC-NAME'])))

    return properties


def create_comments_from_record(record):
    """
    Create Comment elements from compound record data

    Args:
        record: Compound record dictionary

    Returns:
        list: List of Comment objects
    """
    comments = []
    sources = record.get('CREDITS', [])

    source_str = ", ".join(str(s) for s in sources) if isinstance(sources, list) else str(sources)

    if 'COMMENT' in record:
        comment = record['COMMENT']
        if isinstance(comment, list):
            for c in comment:
                comments.append(Comment(value=str(c).strip(), source=source_str))
        elif isinstance(comment, str):
            comments.append(Comment(value=comment.strip(), source=source_str))

    return comments


def create_citation_refs_from_record(record, citation_manager=None):
    """
    Create CitationRef elements from compound record using CitationManager.
    Extracts citations from both CITATIONS field and inline citations in COMMENT field.

    Args:
        record: Compound record dictionary
        citation_manager: CitationManager instance (optional)

    Returns:
        list: List of CitationRef objects (deduplicated)
    """
    if not citation_manager:
        return []

    all_citations = []
    compound_id = record.get('UNIQUE-ID', '')

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
        citation_refs = citation_manager.process_element_citations(compound_id, unique_citations)
        return citation_refs

    return []


def create_enhanced_datanode_from_compound(record, citation_manager=None):
    """
    Convert a compound record to a comprehensive GPML DataNode with citations.

    Args:
        record: Dictionary containing compound data from parsing_utils
        citation_manager: CitationManager instance (optional)

    Returns:
        DataNode: A comprehensive DataNode object
    """
    element_id = record.get('UNIQUE-ID', str(uuid.uuid4()))
    raw_label = record.get('COMMON-NAME', record.get('UNIQUE-ID', 'Unknown Compound'))
    text_label = clean_text_label(raw_label)
    dblinks = record.get('DBLINKS', [])
    db_refs = parse_dblinks(dblinks)

    # Get InChI key from record for primary xref
    inchi_key = record.get('INCHI-KEY', None)

    primary_xref = get_primary_xref(db_refs, inchi_key=inchi_key)

    # Note: Compounds without xref are tracked in statistics, not printed individually
    properties = create_properties_from_record(record)
    comments = create_comments_from_record(record)

    # Get citation references using CitationManager
    citation_refs = create_citation_refs_from_record(record, citation_manager)

    # Use standard metabolite graphics
    graphics = standard_graphics.create_metabolite_graphics(0.0, 0.0)

    datanode = DataNode(
        elementId=element_id,
        textLabel=str(text_label),
        type=DataNodeType.METABOLITE,
        xref=primary_xref,
        states=[],
        graphics=graphics,
        comments=comments,
        properties=properties,
        citationRefs=citation_refs
    )

    return datanode


def create_enhanced_datanodes_from_compounds(compounds_file, citation_manager=None):
    """
    Create comprehensive DataNodes from a compounds file with citation support.

    Args:
        compounds_file: Path to the compounds.dat file
        citation_manager: CitationManager instance (optional)

    Returns:
        list: List of enhanced DataNode objects
    """
    processor = parsing_utils.read_and_parse(compounds_file)

    datanodes = []

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

    # Track compounds without xrefs for examples
    compounds_without_xref = []

    for record in processor.records:
        xref_stats['total_unique_records'] += 1

        # Get primary xref for statistics
        dblinks = record.get('DBLINKS', [])
        db_refs = parse_dblinks(dblinks)
        inchi_key = record.get('INCHI-KEY', None)
        primary_xref = get_primary_xref(db_refs, inchi_key=inchi_key)

        if primary_xref:
            xref_stats['records_with_xref'] += 1
            db_source = primary_xref.dataSource
            xref_stats['xref_by_database'][db_source] = xref_stats['xref_by_database'].get(db_source, 0) + 1
        else:
            # Track compounds without xref (store first 5 examples)
            if len(compounds_without_xref) < 5:
                element_id = record.get('UNIQUE-ID', 'Unknown')
                raw_label = record.get('COMMON-NAME', element_id)
                compounds_without_xref.append(f"{element_id} ({raw_label})")

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
        datanode = create_enhanced_datanode_from_compound(record, citation_manager)
        datanodes.append(datanode)

    # Print statistics
    print("\n" + "="*60)
    print("COMPOUND XREF ANNOTATION STATISTICS")
    print("="*60)
    print(f"Total unique compound records: {xref_stats['total_unique_records']}")
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

    # Show examples of compounds without xrefs
    if compounds_without_xref:
        print(f"\nCompounds without xrefs: {xref_stats['total_unique_records'] - xref_stats['records_with_xref']}")
        print(f"  Examples (first {len(compounds_without_xref)}):")
        for example in compounds_without_xref:
            print(f"    - {example}")

    print("="*60 + "\n")

    return datanodes


def calculate_inchi_key_stats(datanodes):
    """
    Calculate statistics about database annotations for compounds.

    Args:
        datanodes: List of DataNode objects

    Returns:
        dict: Statistics about compound annotations including all databases
    """
    total_compounds = len(datanodes)
    compounds_with_inchi_key = 0
    compounds_with_inchi = 0
    compounds_with_smiles = 0
    compounds_with_primary_inchi_xref = 0

    # Track database annotations from xrefs
    database_counts = {}

    # Track compounds without annotations
    compounds_without_inchi_key = []
    compounds_without_xref = []

    for node in datanodes:
        # Check if primary xref is InChIKey
        has_inchi_key_prop = False
        if node.xref:
            db_source = node.xref.dataSource if node.xref.dataSource else 'Unknown'
            if db_source == 'InChIKey':
                compounds_with_primary_inchi_xref += 1

            # Count this database annotation (only if not empty/unknown)
            if db_source and db_source != 'Unknown':
                if db_source not in database_counts:
                    database_counts[db_source] = 0
                database_counts[db_source] += 1
        else:
            # No xref at all
            compounds_without_xref.append({
                'id': node.elementId,
                'label': node.textLabel
            })

        # Check properties for InChI key and related annotations
        for prop in node.properties:
            if prop.key == 'InChIKey':
                compounds_with_inchi_key += 1
                has_inchi_key_prop = True
            elif prop.key == 'InChI':
                compounds_with_inchi += 1
            elif prop.key == 'SMILES':
                compounds_with_smiles += 1

        # Track compounds without InChI key property
        if not has_inchi_key_prop:
            compounds_without_inchi_key.append({
                'id': node.elementId,
                'label': node.textLabel,
                'xref': f"{node.xref.dataSource}:{node.xref.identifier}" if node.xref else "No xref"
            })

    stats = {
        'total_compounds': total_compounds,
        'compounds_with_inchi_key': compounds_with_inchi_key,
        'compounds_with_inchi': compounds_with_inchi,
        'compounds_with_smiles': compounds_with_smiles,
        'compounds_with_primary_inchi_xref': compounds_with_primary_inchi_xref,
        'inchi_key_percentage': (compounds_with_inchi_key / total_compounds * 100) if total_compounds > 0 else 0,
        'inchi_percentage': (compounds_with_inchi / total_compounds * 100) if total_compounds > 0 else 0,
        'smiles_percentage': (compounds_with_smiles / total_compounds * 100) if total_compounds > 0 else 0,
        'database_annotations': database_counts,
        'compounds_without_inchi_key': compounds_without_inchi_key[:10],  # Keep first 10 examples
        'compounds_without_xref': compounds_without_xref[:10]  # Keep first 10 examples
    }

    return stats


def print_datanode_summary(datanode: DataNode, index):
    """
    Print a summary of a DataNode

    Args:
        datanode: DataNode object
        index: Index number for display
    """
    print(f"\nDataNode {index + 1}:")
    print(f"  ID: {datanode.elementId}")
    print(f"  Label: {datanode.textLabel}")
    print(f"  Type: {datanode.type.value}")

    if datanode.xref:
        print(f"  Primary Xref: {datanode.xref.dataSource}:{datanode.xref.identifier}")

    print(f"  Properties: {len(datanode.properties)} items")
    for prop in datanode.properties[:3]:
        print(f"    - {prop.key}: {prop.value}")
    if len(datanode.properties) > 3:
        print(f"    ... and {len(datanode.properties) - 3} more")

    if datanode.comments:
        print(f"  Comments: {len(datanode.comments)} items")

    if hasattr(datanode, 'citationRefs') and datanode.citationRefs:
        print(f"  Citations: {len(datanode.citationRefs)} references")
        for ref in datanode.citationRefs[:3]:
            print(f"    - {ref.elementRef}")


if __name__ == "__main__":
    # Import CitationManager for testing
    from scripts.build_functions.citation_manager import CitationManager

    # Initialize CitationManager
    citation_manager = CitationManager("pubs.dat")

    # Create datanodes with citation support
    datanodes = create_enhanced_datanodes_from_compounds("compounds.dat", citation_manager)
    print(f"Created {len(datanodes)} enhanced DataNodes")

    # Find compounds with citations for display
    compounds_with_citations = [node for node in datanodes if hasattr(node, 'citationRefs') and node.citationRefs]

    print(f"Found {len(compounds_with_citations)} compounds with citations")

    # Display first few compounds with citations
    for i, node in enumerate(compounds_with_citations[:3]):
        print_datanode_summary(node, i)