from scripts.data_structure.wiki_data_structure import (
    Pathway, Xref, Graphics, Author, Annotation, Citation, Evidence,
    Comment, Property, AnnotationRef, CitationRef, EvidenceRef
)
from scripts.parsing_functions import parsing_utils
from scripts.utils.HTML_cleaner import clean_text_label
from scripts.object2gmpl.gpml_writer import GPMLWriter
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


def parse_pathway_dblinks(dblinks):
    """
    Parse DBLINKS to extract database references for pathways

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


def get_pathway_xref(db_refs):
    """
    Create primary external reference for pathway (typically MetaCyc)

    Args:
        db_refs: Dictionary of database references

    Returns:
        Xref: Primary external reference or None
    """
    db_mapping = {
        'METACYC': 'MetaCyc',
        'KEGG': 'KEGG',
        'REACTOME': 'Reactome',
        'WIKIPATHWAYS': 'WikiPathways',
        'ARACYC': 'AraCyc',
        'PLANTCYC': 'PlantCyc'
    }

    for db_key in ['METACYC', 'KEGG', 'REACTOME', 'WIKIPATHWAYS', 'ARACYC', 'PLANTCYC']:
        if db_key in db_refs:
            return Xref(identifier=db_refs[db_key], dataSource=db_mapping[db_key])

    for db_key, db_id in db_refs.items():
        if db_key in db_mapping:
            return Xref(identifier=db_id, dataSource=db_mapping[db_key])

    return None


def parse_taxonomic_info(record):
    """
    Extract organism and taxonomic information from pathway record.
    Converts ORG-IDs and TAX-IDs to Latin names.

    Args:
        record: Pathway record dictionary

    Returns:
        tuple: (organism_name, taxonomic_ranges)
    """
    organism = None
    taxonomic_ranges = []

    if 'SPECIES' in record:
        species = record['SPECIES']
        if isinstance(species, list):
            # Convert each species ID to Latin name
            converted_species = [convert_to_latin_name(str(s)) for s in species]
            organism = ', '.join(converted_species)
        else:
            organism = convert_to_latin_name(str(species))

    if 'TAXONOMIC-RANGE' in record:
        tax_ranges = record['TAXONOMIC-RANGE']
        if isinstance(tax_ranges, list):
            # Convert each taxonomic range ID to Latin name
            taxonomic_ranges = [convert_to_latin_name(str(t)) for t in tax_ranges]
        else:
            taxonomic_ranges = [convert_to_latin_name(str(tax_ranges))]

    return organism, taxonomic_ranges


def create_pathway_authors_from_record(record):
    """
    Create Author objects from pathway record data

    Args:
        record: Pathway record dictionary

    Returns:
        list: List of Author objects
    """
    authors = []

    if 'CREDITS' in record:
        credits = record['CREDITS']
        if isinstance(credits, list):
            for credit in credits:
                authors.append(Author(name=str(credit)))
        elif isinstance(credits, str):
            authors.append(Author(name=credits))

    return authors


def create_pathway_properties_from_record(record):
    """
    Create Property elements from pathway record data

    Args:
        record: Pathway record dictionary

    Returns:
        list: List of Property objects
    """
    properties = []

    if 'UNIQUE-ID' in record:
        properties.append(Property(key='UniqueID', value=str(record['UNIQUE-ID'])))

    if 'TYPES' in record:
        pathway_types = record['TYPES']
        if isinstance(pathway_types, list):
            for i, ptype in enumerate(pathway_types):
                properties.append(Property(key=f'Type_{i + 1}', value=str(ptype)))
        else:
            properties.append(Property(key='Type', value=str(pathway_types)))

    if 'REACTION-LIST' in record:
        reactions = record['REACTION-LIST']
        if isinstance(reactions, list):
            reaction_count = len(reactions)
            properties.append(Property(key='ReactionCount', value=str(reaction_count)))
            reaction_ids = ', '.join(str(r) for r in reactions)
            properties.append(Property(key='ReactionIDs', value=reaction_ids))
        else:
            properties.append(Property(key='ReactionIDs', value=str(reactions)))

    if 'PREDECESSORS' in record:
        predecessors = record['PREDECESSORS']
        if isinstance(predecessors, list):
            pred_count = len(predecessors)
            properties.append(Property(key='PredecessorCount', value=str(pred_count)))

    if 'PRIMARIES' in record:
        primaries = record['PRIMARIES']
        if isinstance(primaries, list):
            primary_count = len(primaries)
            properties.append(Property(key='PrimaryCompoundCount', value=str(primary_count)))

    if 'PATHWAY-LINKS' in record:
        pathway_links = record['PATHWAY-LINKS']
        if isinstance(pathway_links, list):
            link_text = ', '.join(str(link) for link in pathway_links)
        else:
            link_text = str(pathway_links)
        properties.append(Property(key='PathwayLinks', value=link_text))

    organism, taxonomic_ranges = parse_taxonomic_info(record)
    if organism:
        properties.append(Property(key='Organism', value=organism))

    if taxonomic_ranges:
        for i, tax_range in enumerate(taxonomic_ranges):
            properties.append(Property(key=f'TaxonomicRange_{i + 1}', value=tax_range))

    if 'INSTANCE-NAME-TEMPLATE' in record:
        properties.append(Property(key='InstanceNameTemplate', value=str(record['INSTANCE-NAME-TEMPLATE'])))

    return properties


def create_pathway_comments_from_record(record):
    """
    Create Comment elements from pathway record data

    Args:
        record: Pathway record dictionary

    Returns:
        list: List of Comment objects
    """
    comments = []

    # Add automatic build comment (original BioCyc comments are in description)
    comments.append(Comment(
        value="This GPML file was automatically generated from PlantCyc flatfiles.",
        source="Automated Conversion"
    ))

    return comments


def create_pathway_graphics():
    """
    Create default Graphics object for pathway

    Returns:
        Graphics: Default graphics configuration for pathway
    """
    return Graphics(
        boardWidth=800.0,
        boardHeight=600.0
    )



def create_pathway_citations_from_record_old(record):
    """Original citation creation function - kept for reference."""
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


def create_pathway_citations_from_record(record, citation_manager):
    """
    Create Citation and CitationRef elements from pathway record using CitationManager.
    Extracts citations from both CITATIONS field and inline citations in COMMENT field.

    Args:
        record: Pathway record dictionary
        citation_manager: CitationManager instance

    Returns:
        list: List of CitationRef objects (deduplicated)
    """
    if not citation_manager:
        return []

    all_citations = []
    pathway_id = record.get('UNIQUE-ID', '')

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
        citation_refs = citation_manager.process_element_citations(pathway_id, unique_citations)
        return citation_refs

    return []


def create_enhanced_pathway_from_record(record, citation_manager=None, version=None):
    """
    Convert a pathway record to a comprehensive GPML Pathway object.

    Args:
        record: Dictionary containing pathway data from parsing_utils
        citation_manager: CitationManager instance (optional)
        version: Version string for the pathway (optional)

    Returns:
        Pathway: A comprehensive Pathway object
    """
    element_id = record.get('UNIQUE-ID', str(uuid.uuid4()))
    raw_title = record.get('COMMON-NAME', record.get('UNIQUE-ID', 'Unknown Pathway'))
    title = clean_text_label(raw_title)

    dblinks = record.get('DBLINKS', [])
    db_refs = parse_pathway_dblinks(dblinks)
    xref = get_pathway_xref(db_refs)

    organism, _ = parse_taxonomic_info(record)

    authors = create_pathway_authors_from_record(record)
    properties = create_pathway_properties_from_record(record)
    comments = create_pathway_comments_from_record(record)

    # Use CitationManager if provided, otherwise fall back to old method
    if citation_manager:
        citation_refs = create_pathway_citations_from_record(record, citation_manager)
    else:
        _, citation_refs = create_pathway_citations_from_record_old(record)

    graphics = create_pathway_graphics()

    # Extract description from COMMENT field in record (not from comments list)
    description = None
    if 'COMMENT' in record:
        comment_data = record['COMMENT']
        if isinstance(comment_data, list):
            description = ' '.join(str(c) for c in comment_data)
        else:
            description = str(comment_data)

    pathway = Pathway(
        title=str(title),
        elementId=element_id,
        organism=organism,
        source="PlantCyc",
        version=version,
        license=None,
        xref=xref,
        description=description,
        authors=authors,
        graphics=graphics,
        dataNodes=[],
        interactions=[],
        graphicalLines=[],
        labels=[],
        shapes=[],
        groups=[],
        annotations=[],
        citations=[],  # Will be populated by pathway builder
        evidences=[],
        comments=comments,
        properties=properties,
        annotationRefs=[],
        citationRefs=citation_refs,
        evidenceRefs=[]
    )

    return pathway


def create_enhanced_pathways_from_file(pathways_file, citation_manager=None, version=None):
    """
    Create comprehensive Pathway objects from a pathways file.

    Args:
        pathways_file: Path to the pathways.dat file
        citation_manager: CitationManager instance (optional)
        version: Version string for pathways (optional)

    Returns:
        list: List of enhanced Pathway objects
    """
    processor = parsing_utils.read_and_parse(pathways_file)

    pathways = []

    # Track unique species/taxon
    unique_species = set()
    unique_taxonomic_ranges = set()

    # Citation statistics tracking
    citation_stats = {
        'pathways_with_citations': 0,
        'total_citation_refs': 0,
        'unique_pubmed_ids': set()
    }

    for record in processor.records:
        pathway = create_enhanced_pathway_from_record(record, citation_manager, version)
        pathways.append(pathway)

        # Track species and taxonomic ranges
        organism, taxonomic_ranges = parse_taxonomic_info(record)
        if organism:
            # Split on comma if multiple species are concatenated
            for sp in organism.split(','):
                unique_species.add(sp.strip())
        for tax_range in taxonomic_ranges:
            unique_taxonomic_ranges.add(tax_range)

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

            # Find all inline citations using same pattern as create_pathway_citations_from_record
            inline_pattern = r'\|CITS:\s*\[(\d+)\]\|'
            inline_citations = re.findall(inline_pattern, comment_text)
            all_citations.extend(inline_citations)

        # Deduplicate and track stats
        unique_citations = list(dict.fromkeys(all_citations))

        if unique_citations:
            citation_stats['pathways_with_citations'] += 1
            citation_stats['total_citation_refs'] += len(unique_citations)

            # Extract PubMed IDs
            for citation in unique_citations:
                citation_str = str(citation)
                # PubMed IDs are typically just numbers
                if citation_str.isdigit():
                    citation_stats['unique_pubmed_ids'].add(citation_str)

    # Print statistics
    print("\n" + "="*60)
    print("PATHWAY STATISTICS")
    print("="*60)
    print(f"Total pathways: {len(pathways)}")

    print(f"\nCitation Statistics:")
    print(f"Pathways with citations: {citation_stats['pathways_with_citations']} ({citation_stats['pathways_with_citations']/len(pathways)*100:.1f}%)")
    print(f"Total citation references: {citation_stats['total_citation_refs']}")
    print(f"Unique PubMed IDs: {len(citation_stats['unique_pubmed_ids'])}")

    print(f"\nSpecies/Taxon Information:")
    print(f"Unique species found: {len(unique_species)}")
    print(f"Unique taxonomic ranges found: {len(unique_taxonomic_ranges)}")
    print("="*60 + "\n")

    return pathways


def print_pathway_summary(pathway: Pathway, index):
    """
    Print a summary of a Pathway object

    Args:
        pathway: Pathway object
        index: Index number for display
    """
    print(f"\nPathway {index + 1}:")
    print(f"  ID: {pathway.elementId}")
    print(f"  Title: {pathway.title}")
    print(f"  Organism: {pathway.organism}")
    print(f"  Source: {pathway.source}")

    if pathway.xref:
        print(f"  Primary Xref: {pathway.xref.dataSource}:{pathway.xref.identifier}")
    else:
        print(f"  Primary Xref: None")

    print(f"  Authors ({len(pathway.authors)}):")
    for author in pathway.authors:
        print(f"    - {author.name}")

    print(f"  Properties ({len(pathway.properties)}):")
    for prop in pathway.properties:
        print(f"    - {prop.key}: {prop.value}")

    print(f"  Comments ({len(pathway.comments)}):")
    for comment in pathway.comments:
        print(f"    - {comment.source}: {comment.value}")

    if pathway.citations:
        print(f"  Citations ({len(pathway.citations)}):")
        for citation in pathway.citations:
            xref_str = f"{citation.xref.dataSource}:{citation.xref.identifier}" if citation.xref else "No xref"
            print(f"    - {xref_str}")

    if pathway.description:
        print(f"  Description: {pathway.description}")


if __name__ == "__main__":
    # Import CitationManager only when running as main
    from scripts.build_functions.citation_manager import CitationManager

    # Initialize CitationManager
    citation_manager = CitationManager("pubs.dat")

    # Create pathways with citation support
    pathways = create_enhanced_pathways_from_file("pathways.dat", citation_manager)
    print(f"Created {len(pathways)} enhanced Pathway objects")

    for i, pathway in enumerate(pathways[:3]):
        print_pathway_summary(pathway, i)