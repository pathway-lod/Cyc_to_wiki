from scripts.data_structure.wiki_data_structure import (
    Xref, Annotation, AnnotationType, AnnotationRef
)
import re
from pathlib import Path
import os

# Organism mapping cache
_ORGANISM_MAPPING = None

def load_organism_mapping(mapping_dict=None):
    """Load organism mapping from dictionary or file. Returns dict of dicts."""
    global _ORGANISM_MAPPING

    if mapping_dict is not None:
        _ORGANISM_MAPPING = mapping_dict
        return _ORGANISM_MAPPING

    if _ORGANISM_MAPPING is not None:
        return _ORGANISM_MAPPING

    _ORGANISM_MAPPING = {}

    # Check environment variable first, then fallback to root
    env_path = os.environ.get('ORG_MAPPING_PATH')
    if env_path and Path(env_path).exists():
        mapping_file = Path(env_path)
    else:
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
                    ncbi_id = parts[2] if len(parts) >= 3 else ''
                    
                    _ORGANISM_MAPPING[organism_id] = {
                        'latin_name': latin_name,
                        'ncbi_id': ncbi_id
                    }

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
        info = mapping[organism_id]
        if isinstance(info, dict):
            return info.get('latin_name', organism_id)
        return info

    # Return original if not found
    return organism_id


def get_ncbi_id(organism_id):
    """
    Get NCBI Taxonomy ID for an organism ID.

    Args:
        organism_id: String like 'ORG-5993' or 'TAX-3702'

    Returns:
        str: NCBI ID if found, otherwise None
    """
    if not organism_id:
        return None
        
    mapping = load_organism_mapping()
    organism_id = str(organism_id).strip()
    
    if organism_id in mapping:
        info = mapping[organism_id]
        if isinstance(info, dict):
            return info.get('ncbi_id')
            
    return None


def create_species_annotation(record):
    """
    Create Annotation and AnnotationRef for species/taxonomy.
    Uses NCBI Taxonomy ID for Xref and elementId.

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
            ncbi_id = get_ncbi_id(species_id)
            
            # If we have an NCBI ID, use it for the ID and Xref
            if ncbi_id:
                annotation_id = f"taxonomy_{ncbi_id}"
                xref_id = ncbi_id
                xref_source = "NCBI Taxonomy"
            else:
                # Fallback to internal ID if no NCBI ID available
                annotation_id = f"taxonomy_{species_id}"
                xref_id = species_id
                xref_source = "Taxonomy"
            
            # Sanitize ID
            annotation_id = re.sub(r'[^a-zA-Z0-9_]', '_', annotation_id)
            
            annotation = Annotation(
                elementId=annotation_id,
                value=latin_name,
                type=AnnotationType.TAXONOMY,
                xref=Xref(identifier=xref_id, dataSource=xref_source)
            )
            annotations.append(annotation)
            annotation_refs.append(AnnotationRef(elementRef=annotation_id))

    return annotations, annotation_refs
