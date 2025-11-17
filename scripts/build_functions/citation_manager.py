"""
Citation Manager Module

This module provides comprehensive citation management for BioCyc to GPML conversion,
handling publication loading, citation parsing, and reference management.
"""

import re
from typing import Dict, List, Set, Optional, Tuple
from dataclasses import dataclass, field
from scripts.parsing_functions import parsing_utils
from scripts.data_structure.wiki_data_structure import Citation, CitationRef, Xref


@dataclass
class Publication:
    """Represents a publication with all metadata."""
    unique_id: str
    pubmed_id: Optional[str] = None
    doi_id: Optional[str] = None
    title: Optional[str] = None
    authors: List[str] = field(default_factory=list)
    year: Optional[str] = None
    source: Optional[str] = None
    url: Optional[str] = None
    abstract: Optional[str] = None
    medline_uid: Optional[str] = None
    agricola_id: Optional[str] = None


class CitationManager:
    """
    Manages all citations and publications for BioCyc to GPML conversion.

    This class handles:
    - Loading publications from pubs.dat
    - Parsing various citation formats
    - Managing citation references
    - Creating GPML Citation and CitationRef objects
    """

    def __init__(self, pubs_file: str = "pubs.dat"):
        """
        Initialize the citation manager.

        Args:
            pubs_file: Path to the publications data file
        """
        self.publications: Dict[str, Publication] = {}
        self.pubmed_to_publication: Dict[str, Publication] = {}
        self.citation_aliases: Dict[str, str] = {}  # Maps aliases to unique IDs
        self.citation_objects: Dict[str, Citation] = {}  # Cache of Citation objects
        self.element_citations: Dict[str, Set[str]] = {}  # Maps element IDs to citation IDs

        # Load publications
        self._load_publications(pubs_file)
        self._build_citation_mappings()

    def _load_publications(self, pubs_file: str):
        """
        Load all publications from pubs.dat file.

        Args:
            pubs_file: Path to publications file
        """
        try:
            processor = parsing_utils.read_and_parse(pubs_file)

            for record in processor.records:
                unique_id = record.get('UNIQUE-ID')
                if not unique_id:
                    continue

                # Extract authors list
                authors = []
                if 'AUTHORS' in record:
                    author_data = record['AUTHORS']
                    if isinstance(author_data, list):
                        authors = [str(a) for a in author_data]
                    elif author_data:
                        authors = [str(author_data)]

                # Create Publication object
                pub = Publication(
                    unique_id=unique_id,
                    pubmed_id=str(record.get('PUBMED-ID')) if record.get('PUBMED-ID') else None,
                    doi_id=str(record.get('DOI-ID')) if record.get('DOI-ID') else None,
                    title=str(record.get('TITLE')) if record.get('TITLE') else None,
                    authors=authors,
                    year=str(record.get('YEAR')) if record.get('YEAR') else None,
                    source=str(record.get('SOURCE')) if record.get('SOURCE') else None,
                    url=str(record.get('URL')) if record.get('URL') else None,
                    abstract=str(record.get('ABSTRACT')) if record.get('ABSTRACT') else None,
                    medline_uid=str(record.get('MEDLINE-UID')) if record.get('MEDLINE-UID') else None,
                    agricola_id=str(record.get('AGRICOLA-ID')) if record.get('AGRICOLA-ID') else None
                )

                self.publications[unique_id] = pub

                # Also index by PubMed ID if available
                if pub.pubmed_id:
                    self.pubmed_to_publication[pub.pubmed_id] = pub

            print(f"Loaded {len(self.publications)} publications from {pubs_file}")

        except Exception as e:
            print(f"Warning: Could not load publications from {pubs_file}: {e}")
            print("Citations will be limited to those directly embedded in data files")

    def _build_citation_mappings(self):
        """Build mappings for various citation reference formats."""
        # Build alias mappings for common reference patterns
        for unique_id, pub in self.publications.items():
            # Handle PUB- prefix references
            if unique_id.startswith('PUB-'):
                # Create alias without PUB- prefix
                alias = unique_id[4:].lower()
                self.citation_aliases[alias] = unique_id

                # Also store the full ID as lowercase
                self.citation_aliases[unique_id.lower()] = unique_id

            # Store lowercase version of unique ID
            self.citation_aliases[unique_id.lower()] = unique_id

    def parse_citation_field(self, citation_value: str) -> List[Tuple[str, Optional[str]]]:
        """
        Parse a CITATIONS field value to extract citation ID and evidence code.

        Args:
            citation_value: Raw citation field value

        Returns:
            List of tuples (citation_id, evidence_code)
        """
        citations = []

        if not citation_value:
            return citations

        # Handle list of citations
        if isinstance(citation_value, list):
            for item in citation_value:
                citations.extend(self.parse_citation_field(str(item)))
            return citations

        citation_str = str(citation_value).strip()

        # Parse different citation formats
        # Format 1: PUBMED_ID:EVIDENCE_CODE:OTHER_DATA
        # Format 2: PUBMED_ID
        # Format 3: citation_alias (e.g., "kean81")

        if ':' in citation_str:
            # Split by colon
            parts = citation_str.split(':')
            citation_id = parts[0].strip()
            evidence_code = parts[1].strip() if len(parts) > 1 else None
        else:
            # Simple citation ID
            citation_id = citation_str
            evidence_code = None

        # Clean up citation ID
        citation_id = citation_id.strip()

        if citation_id:
            citations.append((citation_id, evidence_code))

        return citations

    def resolve_citation_id(self, citation_id: str) -> Optional[str]:
        """
        Resolve a citation ID to a publication unique ID.

        Args:
            citation_id: Raw citation identifier

        Returns:
            Resolved publication unique ID or None
        """
        # Check if it's a direct PubMed ID
        if citation_id.isdigit():
            if citation_id in self.pubmed_to_publication:
                return self.pubmed_to_publication[citation_id].unique_id
            # Return as-is for PubMed IDs not in our database
            return f"PUBMED_{citation_id}"

        # Check if it's a known publication unique ID
        if citation_id in self.publications:
            return citation_id

        # Check aliases (case-insensitive)
        citation_lower = citation_id.lower()
        if citation_lower in self.citation_aliases:
            return self.citation_aliases[citation_lower]

        # Try with PUB- prefix
        pub_prefixed = f"PUB-{citation_id.upper()}"
        if pub_prefixed in self.publications:
            return pub_prefixed

        # If not found, create a placeholder ID
        return f"CITATION_{citation_id}"



    def create_citation_object(self, citation_id: str) -> Citation:
        """
        Create a GPML Citation object for a citation ID.

        Args:
            citation_id: Citation identifier

        Returns:
            Citation object for GPML
        """
        # Check cache first
        if citation_id in self.citation_objects:
            return self.citation_objects[citation_id]

        # Resolve the citation ID
        resolved_id = self.resolve_citation_id(citation_id)

        # Create appropriate Xref based on the publication data
        xref = None

        if resolved_id and resolved_id in self.publications:
            pub = self.publications[resolved_id]
            if pub.pubmed_id:
                xref = Xref(identifier=pub.pubmed_id, dataSource="PubMed")
            elif pub.doi_id:
                xref = Xref(identifier=pub.doi_id, dataSource="DOI")
            elif pub.medline_uid:
                xref = Xref(identifier=pub.medline_uid, dataSource="MEDLINE")
            else:
                # Use the unique ID as fallback, but add bibliographic info if available
                # Format: "UNIQUE-ID | Title (Author1, Author2, Year)"
                identifier_parts = [resolved_id]

                if pub.title or pub.authors or pub.year:
                    bib_parts = []
                    if pub.title:
                        # Truncate long titles
                        title = pub.title[:100] + "..." if len(pub.title) > 100 else pub.title
                        bib_parts.append(title)

                    if pub.authors or pub.year:
                        author_year = []
                        if pub.authors:
                            # Take first 3 authors
                            author_list = pub.authors[:3]
                            authors_str = ", ".join(author_list)
                            if len(pub.authors) > 3:
                                authors_str += " et al."
                            author_year.append(authors_str)
                        if pub.year:
                            author_year.append(pub.year)
                        bib_parts.append("(" + ", ".join(author_year) + ")")

                    if bib_parts:
                        identifier_parts.append(" | " + " ".join(bib_parts))

                xref = Xref(identifier="".join(identifier_parts), dataSource="BioCyc")
        elif citation_id.isdigit():
            # Assume it's a PubMed ID
            xref = Xref(identifier=citation_id, dataSource="PubMed")
        else:
            # IMPORTANT: Always create an Xref, even for unknown citations
            # Clean up the citation ID to use as identifier
            clean_id = citation_id.replace('citation_', '').replace('cit_', '')
            xref = Xref(identifier=clean_id, dataSource="PlantCyc")

        # Create Citation object with a sanitized element ID
        element_id = self._sanitize_citation_element_id(resolved_id or citation_id)

        if xref is None:
            # Fallback to ensure we always have an xref
            xref = Xref(identifier=citation_id, dataSource="Unknown")

        citation = Citation(elementId=element_id, xref=xref)

        # Cache the citation object
        self.citation_objects[citation_id] = citation

        return citation

    def _sanitize_citation_element_id(self, citation_id: str) -> str:
        """
        Sanitize citation ID for use as XML element ID.

        Preserves BioCyc publication IDs (PUB-XXXXXX) to maintain references.
        Converts PubMed numeric IDs to valid XML IDs.

        Args:
            citation_id: Raw citation ID (e.g., 'PUB-12695547' or '12695547')

        Returns:
            Sanitized element ID (e.g., 'citation_PUB-12695547' or 'citation_cit_12695547')
        """
        # Keep BioCyc publication IDs intact (PUB-XXXXXX format)
        if citation_id.startswith('PUB-'):
            # Replace hyphen with underscore for XML compatibility
            sanitized = citation_id.replace('-', '_')
            return f"citation_{sanitized}"

        # Handle other prefixes
        if citation_id.startswith('PUBMED_'):
            citation_id = citation_id[7:]
        elif citation_id.startswith('CITATION_'):
            citation_id = citation_id[9:]

        # Replace any invalid characters with underscore
        sanitized = re.sub(r'[^a-zA-Z0-9_]', '_', citation_id)

        # Ensure it doesn't start with a number (for numeric PubMed IDs)
        if sanitized and sanitized[0].isdigit():
            sanitized = f"cit_{sanitized}"

        return f"citation_{sanitized}"

    def process_element_citations(self, element_id: str, citations_field) -> List[CitationRef]:
        """
        Process citations for a specific element and create CitationRef objects.

        Args:
            element_id: ID of the element (pathway, reaction, etc.)
            citations_field: Raw citations field from the record

        Returns:
            List of CitationRef objects (deduplicated)
        """
        citation_refs = []

        if not citations_field:
            return citation_refs

        # Ensure citations_field is a list
        if not isinstance(citations_field, list):
            citations_field = [citations_field]

        # Track citations for this element
        if element_id not in self.element_citations:
            self.element_citations[element_id] = set()

        # Track seen citation element IDs to prevent duplicates
        seen_citation_ids = set()

        for citation_value in citations_field:
            parsed_citations = self.parse_citation_field(citation_value)

            for citation_id, evidence_code in parsed_citations:
                # Create or get citation object
                citation = self.create_citation_object(citation_id)

                # Skip if we've already added this citation
                if citation.elementId in seen_citation_ids:
                    continue

                seen_citation_ids.add(citation.elementId)

                # Track this citation for the element
                self.element_citations[element_id].add(citation.elementId)

                # Create CitationRef
                citation_ref = CitationRef(elementRef=citation.elementId)
                citation_refs.append(citation_ref)

        return citation_refs

    def get_all_citations_for_pathway(self, pathway_elements: List[str]) -> List[Citation]:
        """
        Get all unique Citation objects for a pathway.

        Args:
            pathway_elements: List of element IDs in the pathway

        Returns:
            List of unique Citation objects
        """
        unique_citation_ids = set()

        for element_id in pathway_elements:
            if element_id in self.element_citations:
                unique_citation_ids.update(self.element_citations[element_id])

        # Return Citation objects, ensuring uniqueness by elementId
        citations = []
        seen_element_ids = set()
        for citation_id in unique_citation_ids:
            # Find the citation object by element ID
            for original_id, citation in self.citation_objects.items():
                if citation.elementId == citation_id:
                    # Only add if we haven't seen this elementId before
                    if citation.elementId not in seen_element_ids:
                        citations.append(citation)
                        seen_element_ids.add(citation.elementId)
                    break

        return citations

    def get_publication_info(self, citation_id: str) -> Optional[Publication]:
        """
        Get publication information for a citation.

        Args:
            citation_id: Citation identifier

        Returns:
            Publication object or None
        """
        resolved_id = self.resolve_citation_id(citation_id)
        if resolved_id in self.publications:
            return self.publications[resolved_id]
        return None

    def get_stats(self) -> Dict[str, int]:
        """
        Get statistics about loaded citations.

        Rletseturns:
            Dictionary with citation statistics
        """
        return {
            'total_publications': len(self.publications),
            'pubmed_indexed': len(self.pubmed_to_publication),
            'citation_aliases': len(self.citation_aliases),
            'cached_citations': len(self.citation_objects),
            'elements_with_citations': len(self.element_citations)
        }