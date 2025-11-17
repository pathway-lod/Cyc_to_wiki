"""
Dynamic Property Parser for BioCyc Records
==========================================

Provides utilities for converting BioCyc record fields into GPML Property elements.
"""

from scripts.data_structure.wiki_data_structure import Property


# Field names to skip
SKIP_FIELDS = {
    'COMMON-NAME',      # Already used in textLabel
    'TYPES',            # Internal BioCyc metadata
    'DBLINKS',          # Handled separately for Xrefs
    'CITATIONS',        # Handled separately for Citations
    'COMMENT',          # Handled separately for Comments
    'CREDITS',          # Handled in comment source
}


def normalize_field_name(field_name):
    """
    Convert BioCyc field name to clean property key.

    Examples:
        MOLECULAR-WEIGHT-SEQ -> MolecularWeightSeq
        PI -> Pi
        INCHI-KEY -> InchiKey

    Args:
        field_name: BioCyc field name (e.g., 'MOLECULAR-WEIGHT-SEQ')

    Returns:
        str: Normalized property key
    """
    # Split on hyphens, capitalize each part, join
    parts = field_name.split('-')
    # Capitalize first letter of each part, keep rest as-is for acronyms
    normalized = ''.join(p.capitalize() for p in parts)
    return normalized


def format_value(value):
    """
    Format a field value for storage in a property.

    Args:
        value: Field value (string, list, tuple, etc.)

    Returns:
        str: Formatted value string
    """
    if isinstance(value, list):
        # Handle lists
        if len(value) == 0:
            return None
        elif len(value) == 1:
            return str(value[0]).strip()
        else:
            # Comma-separated list
            return ', '.join(str(v).strip() for v in value if v)

    elif isinstance(value, tuple):
        # Handle tuples (e.g., DBLINKS entries)
        return ', '.join(str(v).strip() for v in value if v)

    elif value is None or str(value).strip() == '':
        return None

    else:
        return str(value).strip()


def create_dynamic_properties(record, skip_fields=None, custom_handlers=None):
    """
    Dynamically parse ALL fields from a BioCyc record into properties.

    This function handles most fields automatically, while allowing custom
    handling for special cases via the custom_handlers dictionary.

    Args:
        record: BioCyc record dictionary
        skip_fields: Additional fields to skip (beyond SKIP_FIELDS)
        custom_handlers: Dict of {field_name: handler_function} for special processing
                        Handler function signature: handler(field_value, record) -> List[Property]

    Returns:
        list: List of Property objects
    """
    properties = []

    # Combine default skip fields with custom ones
    all_skip_fields = SKIP_FIELDS.copy()
    if skip_fields:
        all_skip_fields.update(skip_fields)

    # Initialize custom handlers
    handlers = custom_handlers or {}

    # Process each field in the record
    for field_name, field_value in record.items():
        # Skip if in skip list
        if field_name in all_skip_fields:
            continue

        # Check if there's a custom handler for this field
        if field_name in handlers:
            custom_props = handlers[field_name](field_value, record)
            if custom_props:
                properties.extend(custom_props)
            continue

        # Default handling: convert field name and format value
        formatted_value = format_value(field_value)

        if formatted_value is not None:
            prop_key = normalize_field_name(field_name)
            properties.append(Property(key=prop_key, value=formatted_value))

    return properties


def create_numbered_properties(key_prefix, values):
    """
    Create numbered properties from a list of values.

    Example:
        create_numbered_properties('Synonym', ['ATP', 'adenosine triphosphate'])
        -> [Property(key='Synonym_1', value='ATP'),
            Property(key='Synonym_2', value='adenosine triphosphate')]

    Args:
        key_prefix: Property key prefix
        values: List of values

    Returns:
        list: List of Property objects
    """
    properties = []

    if not isinstance(values, list):
        values = [values]

    for i, value in enumerate(values, 1):
        formatted = format_value(value)
        if formatted:
            properties.append(Property(key=f'{key_prefix}_{i}', value=formatted))

    return properties


# ============================================================================
# Common Custom Handlers
# ============================================================================

def handle_unique_id(field_value, record):
    """Handler for UNIQUE-ID field"""
    return [Property(key='UniqueID', value=str(field_value))]


def handle_synonyms(field_value, record):
    """Handler for SYNONYMS field - creates numbered properties"""
    return create_numbered_properties('Synonym', field_value)


def handle_types(field_value, record):
    """
    Handler for TYPES field - includes it as a property
    """
    formatted = format_value(field_value)
    if formatted:
        return [Property(key='BioCycTypes', value=formatted)]
    return []


def handle_go_terms(field_value, record):
    """Handler for GO-TERMS field"""
    formatted = format_value(field_value)
    if formatted:
        return [Property(key='GOTerms', value=formatted)]
    return []


def handle_chemical_formula(field_value, record):
    """
    Handler for CHEMICAL-FORMULA field - formats tuple structure
    Example: [(H, 2), (O, 1)] -> H 2 O 1
    """
    if isinstance(field_value, list):
        formula = ''
        for part in field_value:
            if isinstance(part, tuple) and len(part) == 2:
                element, count = part
                formula += f"{element} {count} "
            else:
                formula += str(part) + ' '
        return [Property(key='ChemicalFormula', value=formula.strip())]

    formatted = format_value(field_value)
    if formatted:
        return [Property(key='ChemicalFormula', value=formatted)]
    return []
