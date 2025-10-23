from scripts.data_structure.wiki_data_structure import (
    Interaction, Point, Anchor, Graphics, Xref, Property, Comment,
    ArrowHeadType, LineStyle, ConnectorType, AnchorShapeType,
    CitationRef, Citation
)
from scripts.utils.HTML_cleaner import clean_text_label
from scripts.parsing_functions import parsing_utils
import uuid
import re


def parse_regulation_dblinks(dblinks):
    """
    Parse DBLINKS to extract database references for regulations

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


def parse_regulation_mode(mode):
    """
    Parse the MODE field to determine regulation type

    Args:
        mode: MODE field value ('+' for activation, '-' for inhibition)

    Returns:
        tuple: (ArrowHeadType, description string)
    """
    if mode == '+':
        return ArrowHeadType.STIMULATION, "activation"
    elif mode == '-':
        return ArrowHeadType.INHIBITION, "inhibition"
    else:
        # Default to directed arrow for unknown modes
        return ArrowHeadType.DIRECTED, "regulation"


def create_regulation_properties(record):
    """
    Create Property elements from regulation record data

    Args:
        record: Regulation record dictionary

    Returns:
        list: List of Property objects
    """
    properties = []

    # Add ALL database IDs from DBLINKS as individual properties
    dblinks = record.get('DBLINKS', [])
    db_refs = parse_regulation_dblinks(dblinks)
    for db_name, db_id in db_refs.items():
        # Create property key from database name
        prop_key = db_name.replace('-', '_').title().replace('_', '')
        properties.append(Property(key=f'AlternativeId_{prop_key}', value=str(db_id)))

    if 'UNIQUE-ID' in record:
        properties.append(Property(key='RegulationID', value=str(record['UNIQUE-ID'])))

    if 'MODE' in record:
        mode_value = str(record['MODE'])
        mode_description = "activation" if mode_value == '+' else "inhibition" if mode_value == '-' else "unknown"
        properties.append(Property(key='RegulationMode', value=mode_value))
        properties.append(Property(key='RegulationType', value=mode_description))

    if 'MECHANISM' in record:
        properties.append(Property(key='Mechanism', value=str(record['MECHANISM'])))

    if 'PHYSIOLOGICALLY-RELEVANT?' in record:
        properties.append(Property(key='PhysiologicallyRelevant', value=str(record['PHYSIOLOGICALLY-RELEVANT?'])))

    if 'REGULATED-ENTITY' in record:
        properties.append(Property(key='RegulatedEntity', value=str(record['REGULATED-ENTITY'])))

    if 'REGULATOR' in record:
        properties.append(Property(key='Regulator', value=str(record['REGULATOR'])))

    if 'TYPES' in record:
        properties.append(Property(key='RegulationType', value=str(record['TYPES'])))

    return properties


def create_regulation_comments(record):
    """
    Create Comment elements from regulation record data

    Args:
        record: Regulation record dictionary

    Returns:
        list: List of Comment objects
    """
    comments = []

    if 'COMMENT' in record:
        comment_text = record['COMMENT']
        if isinstance(comment_text, list):
            for c in comment_text:
                comments.append(Comment(value=str(c).strip(), source="BioCyc"))
        elif isinstance(comment_text, str):
            comments.append(Comment(value=comment_text.strip(), source="BioCyc"))

    return comments


def create_regulation_citation_refs(record, citation_manager=None):
    """
    Create CitationRef elements from regulation record using CitationManager.
    Extracts citations from both CITATIONS field and inline citations in COMMENT field.

    Args:
        record: Regulation record dictionary
        citation_manager: CitationManager instance (optional)

    Returns:
        list: List of CitationRef objects (deduplicated)
    """
    if not citation_manager:
        return []

    all_citations = []
    regulation_id = record.get('UNIQUE-ID', '')

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
        citation_refs = citation_manager.process_element_citations(regulation_id, unique_citations)
        return citation_refs

    return []


def create_regulation_interaction(record, element_id=None, citation_manager=None):
    """
    Create an Interaction object for a regulation record.

    Regulation interactions represent:
    - Enzyme activity regulation (e.g., allosteric regulation, competitive inhibition)
    - Transcriptional regulation
    - Other regulatory mechanisms

    Args:
        record: Regulation record dictionary
        element_id: Optional element ID
        citation_manager: CitationManager instance (optional)

    Returns:
        dict: Dictionary containing regulation details and connection info
              Returns None if essential fields are missing
    """
    if not element_id:
        element_id = record.get('UNIQUE-ID', str(uuid.uuid4()))

    # Get regulation components
    regulated_entity = record.get('REGULATED-ENTITY')
    regulator = record.get('REGULATOR')
    mode = record.get('MODE')

    # Check if we have essential information
    # If there's no regulated entity or regulator, we can't create a meaningful interaction
    if not regulated_entity and not regulator:
        return None

    # Determine arrow head type based on mode
    arrow_head_type, regulation_type = parse_regulation_mode(mode)

    # Create properties, comments, and citations
    properties = create_regulation_properties(record)
    comments = create_regulation_comments(record)
    citation_refs = create_regulation_citation_refs(record, citation_manager)

    # Create simple xref using the regulation ID
    xref = Xref(identifier=element_id, dataSource='BioCyc')

    # Create waypoints (will be positioned later based on actual nodes)
    waypoints = [
        Point(
            elementId=f"{element_id}_start",
            x=0.0,
            y=0.0,
            arrowHead=ArrowHeadType.UNDIRECTED
        ),
        Point(
            elementId=f"{element_id}_end",
            x=0.0,
            y=0.0,
            arrowHead=arrow_head_type
        )
    ]

    # Create graphics for regulation
    graphics = Graphics(
        lineColor='000000',
        lineStyle=LineStyle.SOLID,
        lineWidth=1.0,
        connectorType=ConnectorType.STRAIGHT
    )

    interaction = Interaction(
        elementId=element_id,
        xref=xref,
        waypoints=waypoints,
        anchors=[],
        graphics=graphics,
        comments=comments,
        properties=properties,
        citationRefs=citation_refs
    )

    return {
        'interaction': interaction,
        'regulated_entity': regulated_entity,
        'regulator': regulator,
        'regulation_type': regulation_type,
        'mode': mode
    }


def create_all_regulation_interactions(regulation_file, citation_manager=None):
    """
    Create all regulation interactions from a regulation file.

    Args:
        regulation_file: Path to the regulation.dat file
        citation_manager: CitationManager instance (optional)

    Returns:
        list: List of regulation interaction dictionaries
    """
    processor = parsing_utils.read_and_parse(regulation_file)

    regulation_interactions = []

    # Statistics tracking
    stats = {
        'total_unique_records': 0,
        'valid_regulations': 0,
        'by_type': {
            'activation': 0,
            'inhibition': 0,
            'unknown': 0
        }
    }

    citation_stats = {
        'records_with_citations': 0,
        'total_citation_refs': 0,
        'unique_pubmed_ids': set()
    }

    # Track examples of unknown regulation types
    unknown_examples = []

    for record in processor.records:
        stats['total_unique_records'] += 1

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

            # Find all inline citations using same pattern as create_regulation_citation_refs
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

        regulation_data = create_regulation_interaction(record, citation_manager=citation_manager)
        if regulation_data:  # Only add if we got valid data
            regulation_interactions.append(regulation_data)
            stats['valid_regulations'] += 1

            # Track by regulation type
            reg_type = regulation_data.get('regulation_type', 'unknown')
            if reg_type in stats['by_type']:
                stats['by_type'][reg_type] += 1
            else:
                stats['by_type']['unknown'] += 1

            # Track examples of unknown regulation types
            if reg_type == 'unknown' and len(unknown_examples) < 5:
                unique_id = record.get('UNIQUE-ID', 'N/A')
                mode = record.get('MODE', 'N/A')
                regulator = regulation_data.get('regulator', 'N/A')
                regulated = regulation_data.get('regulated_entity', 'N/A')
                unknown_examples.append(f"{unique_id} (MODE: {mode}, {regulator} -> {regulated})")

    # Print statistics
    print("\n" + "="*60)
    print("REGULATION INTERACTION STATISTICS")
    print("="*60)
    print(f"Total unique regulation records: {stats['total_unique_records']}")
    print(f"Valid regulations created: {stats['valid_regulations']}")
    if stats['total_unique_records'] > 0:
        print(f"Valid rate: {stats['valid_regulations']/stats['total_unique_records']*100:.1f}%")
    print(f"\nRegulations by type:")
    for reg_type in sorted(stats['by_type'].keys()):
        count = stats['by_type'][reg_type]
        percentage = (count / stats['valid_regulations'] * 100) if stats['valid_regulations'] > 0 else 0
        print(f"  {reg_type}: {count} ({percentage:.1f}%)")

    # Print examples of unknown regulation types
    if unknown_examples:
        print(f"\nExamples of unknown regulation types (first {len(unknown_examples)}):")
        for example in unknown_examples:
            print(f"    - {example}")

    print(f"\nCitation Statistics:")
    print(f"Records with citations: {citation_stats['records_with_citations']} ({citation_stats['records_with_citations']/stats['total_unique_records']*100:.1f}%)")
    print(f"Total citation references: {citation_stats['total_citation_refs']}")
    print(f"Unique PubMed IDs: {len(citation_stats['unique_pubmed_ids'])}")
    print("="*60 + "\n")

    return regulation_interactions


def print_regulation_interaction_summary(regulation_data, index):
    """
    Print a summary of a regulation interaction

    Args:
        regulation_data: Dictionary containing regulation interaction data
        index: Index number for display
    """
    interaction = regulation_data['interaction']

    print(f"\nRegulation Interaction {index + 1}:")
    print(f"  ID: {interaction.elementId}")

    if interaction.xref:
        print(f"  Primary Xref: {interaction.xref.dataSource}:{interaction.xref.identifier}")

    print(f"  Regulation Type: {regulation_data['regulation_type']}")
    print(f"  Mode: {regulation_data['mode']}")
    print(f"  Regulated Entity: {regulation_data['regulated_entity']}")
    print(f"  Regulator: {regulation_data['regulator']}")

    print(f"  Properties: {len(interaction.properties)} items")
    for prop in interaction.properties[:3]:
        value_preview = prop.value[:50] + '...' if len(prop.value) > 50 else prop.value
        print(f"    - {prop.key}: {value_preview}")


if __name__ == "__main__":
    regulations = create_all_regulation_interactions("regulation.dat")
    print(f"Created {len(regulations)} regulation interactions")

    for i, regulation_data in enumerate(regulations[:5]):
        print_regulation_interaction_summary(regulation_data, i)
