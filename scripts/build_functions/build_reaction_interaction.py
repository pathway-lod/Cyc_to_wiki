from scripts.data_structure.wiki_data_structure import (
    Interaction, Point, Anchor, Graphics, Xref, Property, Comment,
    ArrowHeadType, LineStyle, ConnectorType, AnchorShapeType,
    CitationRef, Citation
)
from scripts.utils.HTML_cleaner import clean_text_label
from scripts.utils.property_parser import create_dynamic_properties
from scripts.parsing_functions import parsing_utils
import uuid
import re


def parse_reaction_dblinks(dblinks):
    """
    Parse DBLINKS to extract database references for reactions

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


def get_primary_reaction_xref(db_refs, unique_id=None, ec_number=None):
    """
    Get the primary external reference for a reaction using Rhea as priority

    Args:
        db_refs: Dictionary of database references
        unique_id: The UNIQUE-ID from the reaction record
        ec_number: EC number if available

    Returns:
        Xref: Primary external reference or None
    """
    # Priority 1: Rhea from DBLINKS
    if 'RHEA' in db_refs:
        return Xref(identifier=db_refs['RHEA'], dataSource='Rhea')

    # Priority 2: EC number
    if ec_number:
        return Xref(identifier=str(ec_number), dataSource='Enzyme Nomenclature')

    # Priority 3: Other database references
    db_mapping = {
        'METACYC': 'MetaCyc',
        'KEGG': 'KEGG Reaction',
        'BIGG': 'BiGG',
        'SEED': 'ModelSEED',
        'METANETX': 'MetaNetX'
    }
    for db_key in ['METACYC', 'KEGG', 'BIGG', 'SEED', 'METANETX']:
        if db_key in db_refs:
            return Xref(identifier=db_refs[db_key], dataSource=db_mapping[db_key])

    # Priority 4: Use UNIQUE-ID as final fallback
    if unique_id:
        # Database mapping
        source_db_mapping = {
            'METACYC': 'MetaCyc',
            'ECOCYC': 'EcoCyc',
            'XMETDB': 'XMetDB',
            'BIOCYC': 'BioCyc',
            'PLANTCYC' : 'PlantCyc'
        }

        # Determine database source from DBLINKS
        source_db = None
        for db_key in ['METACYC', 'ECOCYC', 'BIOCYC', 'XMETDB']:
            if db_key in db_refs:
                source_db = db_key
                break

        # If no specific source found, default to MetaCyc
        if not source_db:
            source_db = 'PLANTCYC'

        datasource = source_db_mapping.get(source_db, 'PlantCyc')
        return Xref(identifier=unique_id, dataSource=datasource)

    return None


def create_reaction_properties(record):
    """
    Create Property elements from reaction record data

    Args:
        record: Reaction record dictionary

    Returns:
        list: List of Property objects
    """
    properties = []

    # Add ALL database IDs from DBLINKS as individual properties
    dblinks = record.get('DBLINKS', [])
    db_refs = parse_reaction_dblinks(dblinks)
    for db_name, db_id in db_refs.items():
        # Create property key from database name
        prop_key = db_name.replace('-', '_').title().replace('_', '')
        properties.append(Property(key=f'AlternativeId_{prop_key}', value=str(db_id)))

    if 'UNIQUE-ID' in record:
        properties.append(Property(key='UniqueID', value=str(record['UNIQUE-ID'])))

    if 'COMMON-NAME' in record:
        properties.append(Property(key='ReactionName', value=str(record['COMMON-NAME'])))

    if 'EC-NUMBER' in record:
        ec_numbers = record['EC-NUMBER']
        if isinstance(ec_numbers, list):
            for i, ec in enumerate(ec_numbers):
                properties.append(Property(key=f'EC_Number_{i + 1}', value=str(ec)))
        else:
            properties.append(Property(key='EC_Number', value=str(ec_numbers)))

    if 'OFFICIAL-EC?' in record:
        properties.append(Property(key='OfficialEC', value=str(record['OFFICIAL-EC?'])))

    if 'DELTAG0' in record:
        properties.append(Property(key='StandardGibbsFreeEnergy', value=str(record['DELTAG0'])))

    if 'SPONTANEOUS?' in record:
        properties.append(Property(key='Spontaneous', value=str(record['SPONTANEOUS?'])))

    if 'ORPHAN?' in record:
        properties.append(Property(key='Orphan', value=str(record['ORPHAN?'])))

    if 'IN-PATHWAY' in record:
        pathways = record['IN-PATHWAY']
        if isinstance(pathways, list):
            pathway_list = ', '.join(str(p) for p in pathways)
            properties.append(Property(key='Pathways', value=pathway_list))
            properties.append(Property(key='PathwayCount', value=str(len(pathways))))
        else:
            properties.append(Property(key='Pathways', value=str(pathways)))

    #build reaction equation
    equation = build_reaction_equation(record)
    if equation:
        properties.append(Property(key='ReactionEquation', value=equation))

    if 'SIGNAL' in record:
        signal = record['SIGNAL']
        if isinstance(signal, list):
            signal_text = ', '.join(str(s) for s in signal)
        else:
            signal_text = str(signal)
        properties.append(Property(key='Signal', value=signal_text))

    if 'SPECIES' in record:
        species = record['SPECIES']
        if isinstance(species, list):
            species_text = ', '.join(str(s) for s in species)
        else:
            species_text = str(species)
        properties.append(Property(key='Species', value=species_text))

    if 'REACTION-DIRECTION' in record:
        direction = str(record['REACTION-DIRECTION'])
        properties.append(Property(key='ReactionDirection', value=direction))

    if 'EQUILIBRIUM-CONSTANT' in record:
        properties.append(Property(key='EquilibriumConstant', value=str(record['EQUILIBRIUM-CONSTANT'])))

    if 'PHI-MINUS-PHI-ZERO' in record:
        properties.append(Property(key='PhiMinusPhiZero', value=str(record['PHI-MINUS-PHI-ZERO'])))

    if 'STD-REDUCTION-POTENTIAL' in record:
        properties.append(Property(key='StandardReductionPotential', value=str(record['STD-REDUCTION-POTENTIAL'])))

    if 'REACTION-BALANCE-STATUS' in record:
        properties.append(Property(key='BalanceStatus', value=str(record['REACTION-BALANCE-STATUS'])))

    if 'PHYSIOLOGICALLY-RELEVANT?' in record:
        properties.append(Property(key='PhysiologicallyRelevant', value=str(record['PHYSIOLOGICALLY-RELEVANT?'])))

    if 'ENZYMATIC-REACTION' in record:
        enzrxns = record['ENZYMATIC-REACTION']
        if isinstance(enzrxns, list):
            enzrxn_text = ', '.join(str(e) for e in enzrxns)
            properties.append(Property(key='EnzymaticReactions', value=enzrxn_text))
            properties.append(Property(key='EnzymaticReactionCount', value=str(len(enzrxns))))
        else:
            properties.append(Property(key='EnzymaticReaction', value=str(enzrxns)))

    if 'SYNONYMS' in record:
        synonyms = record['SYNONYMS']
        if isinstance(synonyms, list):
            for i, synonym in enumerate(synonyms):
                properties.append(Property(key=f'Synonym_{i + 1}', value=str(synonym).strip()))
        elif isinstance(synonyms, str):
            properties.append(Property(key='Synonym_1', value=synonyms.strip()))

    if 'SYSTEMATIC-NAME' in record:
        properties.append(Property(key='SystematicName', value=str(record['SYSTEMATIC-NAME'])))

    if 'TYPES' in record:
        types = record['TYPES']
        if isinstance(types, list):
            types_text = ', '.join(str(t) for t in types)
            properties.append(Property(key='ReactionTypes', value=types_text))
        else:
            properties.append(Property(key='ReactionType', value=str(types)))

    # Dynamic parsing for all other fields
    skip_fields = {
        'UNIQUE-ID', 'COMMON-NAME', 'DBLINKS', 'EC-NUMBER', 'OFFICIAL-EC?',
        'DELTAG0', 'SPONTANEOUS?', 'ORPHAN?', 'IN-PATHWAY', 'REACTION-DIRECTION',
        'REACTION-BALANCE-STATUS', 'LEFT', 'RIGHT', 'PHYSIOLOGICALLY-RELEVANT?',
        'STD-REDUCTION-POTENTIAL', 'ENZYME-REACTION-CLASS', 'ENZYMATIC-REACTION',
        'CATALYZED-BY', 'RXN-LOCATIONS', 'TYPES',
        'CITATIONS', 'COMMENT', 'CREDITS',
    }
    dynamic_props = create_dynamic_properties(record, skip_fields=skip_fields)
    properties.extend(dynamic_props)

    return properties


def create_reaction_comments(record):
    """
    Create Comment elements from reaction record data

    Args:
        record: Reaction record dictionary

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


def create_reaction_citation_refs(record, citation_manager=None):
    """
    Create CitationRef elements from reaction record using CitationManager.
    Extracts citations from both CITATIONS field and inline citations in COMMENT field.

    Args:
        record: Reaction record dictionary
        citation_manager: CitationManager instance (optional)

    Returns:
        list: List of CitationRef objects (deduplicated)
    """
    if not citation_manager:
        return []

    all_citations = []
    reaction_id = record.get('UNIQUE-ID', '')

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
        citation_refs = citation_manager.process_element_citations(reaction_id, unique_citations)
        return citation_refs

    return []

def parse_reaction_components(record):
    """
    Parse LEFT and RIGHT components of a reaction

    Args:
        record: Reaction record dictionary

    Returns:
        tuple: (reactants_list, products_list) where each item is a dict with compound_id and coefficient
    """
    reactants = []
    products = []

    if 'LEFT' in record:
        left = record['LEFT']
        if isinstance(left, list):
            for reactant in left:
                if isinstance(reactant, tuple) and len(reactant) >= 2:
                    coefficient, compound_id = reactant[0], reactant[1]
                    reactants.append({
                        'compound_id': str(compound_id),
                        'coefficient': float(coefficient) if coefficient else 1.0
                    })
                else:
                    reactants.append({
                        'compound_id': str(reactant),
                        'coefficient': 1.0
                    })
        else:
            reactants.append({
                'compound_id': str(left),
                'coefficient': 1.0
            })

    if 'RIGHT' in record:
        right = record['RIGHT']
        if isinstance(right, list):
            for product in right:
                if isinstance(product, tuple) and len(product) >= 2:
                    coefficient, compound_id = product[0], product[1]
                    products.append({
                        'compound_id': str(compound_id),
                        'coefficient': float(coefficient) if coefficient else 1.0
                    })
                else:
                    products.append({
                        'compound_id': str(product),
                        'coefficient': 1.0
                    })
        else:
            products.append({
                'compound_id': str(right),
                'coefficient': 1.0
            })

    return reactants, products


def build_reaction_equation(record):
    """
    Build a human-readable reaction equation

    Args:
        record: Reaction record dictionary

    Returns:
        str: Formatted reaction equation
    """
    reactants, products = parse_reaction_components(record)

    # Format reactants
    reactant_parts = []
    for r in reactants:
        compound_label = clean_text_label(r['compound_id'])
        if r['coefficient'] != 1.0:
            reactant_parts.append(f"{r['coefficient']} {compound_label}")
        else:
            reactant_parts.append(compound_label)

    # Format products
    product_parts = []
    for p in products:
        compound_label = clean_text_label(p['compound_id'])
        if p['coefficient'] != 1.0:
            product_parts.append(f"{p['coefficient']} {compound_label}")
        else:
            product_parts.append(compound_label)

    # Determine arrow type
    is_reversible = record.get('REACTION-DIRECTION') == 'REVERSIBLE'
    arrow = " <=> " if is_reversible else " => "

    return " + ".join(reactant_parts) + arrow + " + ".join(product_parts)


def determine_arrow_head_type(record):
    """
    Determine the appropriate arrow head type for a reaction

    Args:
        record: Reaction record dictionary

    Returns:
        ArrowHeadType: The appropriate arrow head type
    """
    # Check if reaction is reversible
    if record.get('REACTION-DIRECTION') == 'REVERSIBLE':
        return ArrowHeadType.CONVERSION  # Double arrow if reversible

    if 'ENZYMATIC-REACTION' in record:
        return ArrowHeadType.CATALYSIS

    # use directed arrow for non-reversible reactions
    return ArrowHeadType.DIRECTED


def create_reaction_interaction_with_anchors(record, element_id=None, citation_manager=None):
    """
    Create an Interaction with anchor points for connecting multiple reactants/products.

    Args:
        record: Reaction record dictionary
        element_id: Optional element ID
        citation_manager: CitationManager instance (optional)

    Returns:
        dict: Dictionary containing interaction details and connection info
    """
    if not element_id:
        element_id = record.get('UNIQUE-ID', str(uuid.uuid4()))

    reactants, products = parse_reaction_components(record)

    # Create waypoints
    waypoints = [
        Point(
            elementId=f"{element_id}_wp1",
            x=200.0,
            y=300.0,
            arrowHead=ArrowHeadType.UNDIRECTED
        ),
        Point(
            elementId=f"{element_id}_wp2",
            x=400.0,
            y=300.0,
            arrowHead=determine_arrow_head_type(record)
        )
    ]

    # Create anchors for multiple reactants/products
    anchors = []
    if len(reactants) > 1:
        for i in range(len(reactants)):
            anchor = Anchor(
                elementId=f"{element_id}_reactant_anchor_{i}",
                position=0.25,
                shapeType=AnchorShapeType.CIRCLE
            )
            anchors.append(anchor)

    if len(products) > 1:
        for i in range(len(products)):
            anchor = Anchor(
                elementId=f"{element_id}_product_anchor_{i}",
                position=0.75,
                shapeType=AnchorShapeType.CIRCLE
            )
            anchors.append(anchor)

    # Parse database links
    dblinks = record.get('DBLINKS', [])
    db_refs = parse_reaction_dblinks(dblinks)
    unique_id = record.get('UNIQUE-ID')
    ec_number = record.get('EC-NUMBER')
    if ec_number and isinstance(ec_number, list):
        ec_number = ec_number[0]

    xref = get_primary_reaction_xref(db_refs, unique_id, ec_number)
    properties = create_reaction_properties(record)
    comments = create_reaction_comments(record)

    # Get citation references using CitationManager
    citation_refs = create_reaction_citation_refs(record, citation_manager)

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
        anchors=anchors,
        graphics=graphics,
        comments=comments,
        properties=properties,
        citationRefs=citation_refs
    )

    return {
        'interaction': interaction,
        'reactants': reactants,
        'products': products,
        'reactant_anchors': [a.elementId for a in anchors if 'reactant' in a.elementId],
        'product_anchors': [a.elementId for a in anchors if 'product' in a.elementId]
    }


def create_all_reaction_interactions(reactions_file, citation_manager=None):
    """
    Create all reaction interactions from a reactions file.

    Args:
        reactions_file: Path to the reactions.dat file
        citation_manager: CitationManager instance (optional)

    Returns:
        list: List of reaction interaction dictionaries
    """
    processor = parsing_utils.read_and_parse(reactions_file)

    reaction_interactions = []

    # Statistics tracking
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

    for record in processor.records:
        xref_stats['total_unique_records'] += 1

        # Get primary xref for statistics
        dblinks = record.get('DBLINKS', [])
        db_refs = parse_reaction_dblinks(dblinks)
        unique_id = record.get('UNIQUE-ID')
        ec_number = record.get('EC-NUMBER')
        if ec_number and isinstance(ec_number, list):
            ec_number = ec_number[0]
        primary_xref = get_primary_reaction_xref(db_refs, unique_id, ec_number)

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

            # Find all inline citations using same pattern as create_reaction_citation_refs
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

        # Create reaction interaction
        reaction_data = create_reaction_interaction_with_anchors(record,
                                                                 citation_manager=citation_manager)
        reaction_interactions.append(reaction_data)

    # Print statistics
    print("\n" + "="*60)
    print("REACTION XREF ANNOTATION STATISTICS")
    print("="*60)
    print(f"Total unique reaction records: {xref_stats['total_unique_records']}")
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

    return reaction_interactions


def print_reaction_interaction_summary(reaction_data, index):
    """
    Print a summary of a reaction interaction

    Args:
        reaction_data: Dictionary containing reaction interaction data
        index: Index number for display
    """
    interaction = reaction_data['interaction']

    print(f"\nReaction Interaction {index + 1}:")
    print(f"  ID: {interaction.elementId}")

    if interaction.xref:
        print(f"  Primary Xref: {interaction.xref.dataSource}:{interaction.xref.identifier}")

    print(f"  Properties: {len(interaction.properties)} items")
    for prop in interaction.properties[:3]:
        value_preview = prop.value[:50] + '...' if len(prop.value) > 50 else prop.value
        print(f"    - {prop.key}: {value_preview}")

    print(f"  Reactants: {len(reaction_data['reactants'])}")
    for r in reaction_data['reactants'][:3]:
        print(f"    - {r['coefficient']} x {r['compound_id']}")

    print(f"  Products: {len(reaction_data['products'])}")
    for p in reaction_data['products'][:3]:
        print(f"    - {p['coefficient']} x {p['compound_id']}")

    print(f"  Anchors: {len(interaction.anchors)} total")
    print(f"    - Reactant anchors: {len(reaction_data['reactant_anchors'])}")
    print(f"    - Product anchors: {len(reaction_data['product_anchors'])}")


if __name__ == "__main__":
    reactions = create_all_reaction_interactions("reactions.dat")
    print(f"Created {len(reactions)} reaction interactions")

    for i, reaction_data in enumerate(reactions[:3]):
        print_reaction_interaction_summary(reaction_data, i)