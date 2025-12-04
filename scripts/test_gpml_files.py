#!/usr/bin/env python3
"""
GPML File Validation Script

Validates GPML (Graphical Pathway Markup Language) files for:
1. XML schema compliance (GPML 2021)
2. Duplicate element IDs (including specific checks for duplicate interaction IDs)
3. Unsanitized IDs (invalid XML characters)
4. Missing references to compounds in interactions
5. Missing references to groups and citation elements
6. All element references and cross-references
7. Extra checks based on experiences with Pathvisio crashes

Usage:
    python test_gpml_files.py <gpml_file_or_directory>

    # Single file
    python test_gpml_files.py pathway.gpml

    # Directory
    python test_gpml_files.py biocyc_pathways_*/
"""

import os
import sys
import json
import re
from pathlib import Path
from typing import Dict, List, Set, Tuple, Any
from dataclasses import dataclass, field
from xml.etree import ElementTree as ET
from collections import defaultdict


# GPML Namespace
GPML_NS = "http://pathvisio.org/GPML/2021"
NS = {"gpml": GPML_NS}

# ID validation regex - must start with letter or underscore, contain only alphanumeric and underscore
# !Still don't know if this is intendent behaviour or bug in Pathvisio!!
VALID_ID_PATTERN = re.compile(r'^[a-zA-Z_][a-zA-Z0-9_]*$')

# Arrow head types
VALID_ARROW_HEADS = {
    "Undirected", "Directed", "Conversion", "Inhibition",
    "Catalysis", "Stimulation", "Binding", "Translocation",
    "TranscriptionTranslation"
}

# Valid data node types
VALID_DATANODE_TYPES = {
    "GeneProduct", "Protein", "Metabolite", "Complex",
    "DNA", "RNA", "PATHWAY", "DISEASE", "PHENOTYPE", "ALIAS", "EVENT"
}

# Valid group types
VALID_GROUP_TYPES = {
    "Group", "Complex", "Pathway", "Transparent", "Analog", "Paralog"
}


@dataclass
class ValidationError:
    """Represents a validation error"""
    file_path: str
    error_type: str
    message: str
    element_id: str = ""
    line_number: int = 0


@dataclass
class ValidationReport:
    """Contains all validation results"""
    file_path: str
    xml_valid: bool = False
    errors: List[ValidationError] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    stats: Dict[str, Any] = field(default_factory=dict)

    def add_error(self, error_type: str, message: str, element_id: str = ""):
        """Add an error to the report"""
        self.errors.append(ValidationError(
            file_path=self.file_path,
            error_type=error_type,
            message=message,
            element_id=element_id
        ))

    def add_warning(self, message: str):
        """Add a warning to the report"""
        self.warnings.append(message)

    def is_valid(self) -> bool:
        """Return True if no errors found"""
        return len(self.errors) == 0


class GPMLValidator:
    """Validates GPML files"""

    def __init__(self):
        self.report = None
        self.all_element_ids = set()
        self.datanode_ids = set()
        self.interaction_ids = set()
        self.group_ids = set()
        self.citation_ids = set()
        self.state_ids = set()
        self.point_ids = set()
        self.anchor_ids = set()

    def validate_file(self, file_path: str) -> ValidationReport:
        """Validate a single GPML file"""
        self.report = ValidationReport(file_path=file_path)
        self.all_element_ids = set()
        self.datanode_ids = set()
        self.interaction_ids = set()
        self.group_ids = set()
        self.citation_ids = set()
        self.state_ids = set()
        self.point_ids = set()
        self.anchor_ids = set()

        # Try to parse XML
        try:
            tree = ET.parse(file_path)
            root = tree.getroot()
            self.report.xml_valid = True
        except ET.ParseError as e:
            self.report.add_error("XML_PARSE_ERROR", f"XML parsing failed: {str(e)}")
            return self.report
        except FileNotFoundError:
            self.report.add_error("FILE_NOT_FOUND", f"File not found: {file_path}")
            return self.report
        except Exception as e:
            self.report.add_error("XML_ERROR", f"Error reading file: {str(e)}")
            return self.report

        # Validate root element
        if root.tag != f"{{{GPML_NS}}}Pathway":
            self.report.add_error("INVALID_ROOT", f"Root element is not 'Pathway', got '{root.tag}'")

        # Collect statistics
        self.report.stats = {
            "datanodes": len(root.findall("gpml:DataNodes/gpml:DataNode", NS)),
            "interactions": len(root.findall("gpml:Interactions/gpml:Interaction", NS)),
            "groups": len(root.findall("gpml:Groups/gpml:Group", NS)),
            "citations": len(root.findall("gpml:Citations/gpml:Citation", NS)),
        }

        # Run validation checks
        self._validate_schema(root)
        self._validate_organism(root)
        self._validate_duplicate_ids(root)
        self._validate_id_format(root)
        self._validate_references(root)
        self._validate_required_fields(root)
        self._validate_graphics(root)
        self._validate_xrefs(root)
        self._validate_circular_references(root)
        self._validate_missing_elements(root)

        return self.report

    def _validate_schema(self, root: ET.Element):
        """Validate XML schema compliance"""
        # Check required pathway attributes
        if root.get("title") is None:
            self.report.add_warning("Missing 'title' attribute on Pathway")

        graphics = root.find("gpml:Graphics", NS)
        if graphics is None:
            self.report.add_error("SCHEMA_ERROR", "Missing Graphics element in Pathway")
        else:
            if graphics.get("boardWidth") is None:
                self.report.add_error("SCHEMA_ERROR", "Missing 'boardWidth' in Pathway Graphics")
            if graphics.get("boardHeight") is None:
                self.report.add_error("SCHEMA_ERROR", "Missing 'boardHeight' in Pathway Graphics")

    def _validate_organism(self, root: ET.Element):
        """Validate organism attribute on Pathway element"""
        organism = root.get("organism")

        # Check if organism attribute exists
        if organism is None or organism.strip() == "":
            self.report.add_error(
                "MISSING_ORGANISM",
                "Pathway element missing 'organism' attribute or organism is empty"
            )
            return

        # Check for multiple organisms
        if "," in organism:
            organisms = [o.strip() for o in organism.split(",")]
            self.report.add_error(
                "MULTIPLE_ORGANISMS",
                f"Pathway has multiple organisms ({len(organisms)}): {organism}. "
                f"Each GPML file must have exactly one organism."
            )

    def _validate_duplicate_ids(self, root: ET.Element):
        """Check for duplicate element IDs across all elements"""
        id_locations = defaultdict(list)

        # Track interaction IDs separately to check for duplicates within interactions
        interaction_id_list = []

        # Collect all elements with IDs
        for element in root.iter():
            elem_id = element.get("elementId")
            if elem_id:
                tag_name = element.tag.split("}")[-1] if "}" in element.tag else element.tag
                id_locations[elem_id].append(tag_name)
                self.all_element_ids.add(elem_id)

                # Track interaction IDs in order
                if tag_name == "Interaction":
                    interaction_id_list.append(elem_id)

        # Find duplicates across all elements
        for elem_id, tag_names in id_locations.items():
            if len(tag_names) > 1:
                self.report.add_error(
                    "DUPLICATE_ID",
                    f"Duplicate elementId: '{elem_id}' found in {tag_names}",
                    elem_id
                )

            # Track by element type
            if tag_names and tag_names[0] == "DataNode":
                self.datanode_ids.add(elem_id)
            elif tag_names and tag_names[0] == "Interaction":
                self.interaction_ids.add(elem_id)
            elif tag_names and tag_names[0] == "Group":
                self.group_ids.add(elem_id)
            elif tag_names and tag_names[0] == "Citation":
                self.citation_ids.add(elem_id)
            elif tag_names and tag_names[0] == "State":
                self.state_ids.add(elem_id)
            elif tag_names and tag_names[0] == "Point":
                self.point_ids.add(elem_id)
            elif tag_names and tag_names[0] == "Anchor":
                self.anchor_ids.add(elem_id)

        # Check specifically for duplicate interaction IDs
        seen_interaction_ids = set()
        for interaction_id in interaction_id_list:
            if interaction_id in seen_interaction_ids:
                self.report.add_error(
                    "DUPLICATE_INTERACTION_ID",
                    f"Duplicate Interaction elementId: '{interaction_id}'",
                    interaction_id
                )
            seen_interaction_ids.add(interaction_id)

    def _validate_id_format(self, root: ET.Element):
        """Check for unsanitized IDs (invalid XML characters)"""
        for element in root.iter():
            elem_id = element.get("elementId")
            if elem_id:
                # Check for spaces
                if " " in elem_id:
                    self.report.add_error(
                        "INVALID_ID_FORMAT",
                        f"Element ID contains spaces: '{elem_id}'",
                        elem_id
                    )

                # Check for problematic XML special characters
                special_chars = set(elem_id) & set("<>\"'&/\\")
                if special_chars:
                    self.report.add_error(
                        "INVALID_ID_FORMAT",
                        f"Element ID contains special characters {special_chars}: '{elem_id}'",
                        elem_id
                    )

    def _validate_references(self, root: ET.Element):
        """Validate all element references (groupRef, elementRef, etc.)"""

        # Validate DataNode groupRef references
        for datanode in root.findall("gpml:DataNodes/gpml:DataNode", NS):
            datanode_id = datanode.get("elementId", "UNKNOWN")
            group_ref = datanode.get("groupRef")
            if group_ref and group_ref not in self.group_ids and group_ref not in self.all_element_ids:
                self.report.add_error(
                    "MISSING_GROUP_REF",
                    f"DataNode '{datanode_id}' references missing group '{group_ref}'",
                    datanode_id
                )

            # Validate CitationRef references
            for citation_ref in datanode.findall("gpml:CitationRef", NS):
                citation_id = citation_ref.get("elementRef")
                if citation_id and citation_id not in self.citation_ids and citation_id not in self.all_element_ids:
                    self.report.add_error(
                        "MISSING_CITATION_REF",
                        f"DataNode '{datanode_id}' references missing citation '{citation_id}'",
                        datanode_id
                    )

        # Validate Interaction groupRef and Point references
        for interaction in root.findall("gpml:Interactions/gpml:Interaction", NS):
            interaction_id = interaction.get("elementId", "UNKNOWN")
            group_ref = interaction.get("groupRef")
            if group_ref and group_ref not in self.group_ids and group_ref not in self.all_element_ids:
                self.report.add_error(
                    "MISSING_GROUP_REF",
                    f"Interaction '{interaction_id}' references missing group '{group_ref}'",
                    interaction_id
                )

            # Validate CitationRef references
            for citation_ref in interaction.findall("gpml:CitationRef", NS):
                citation_id = citation_ref.get("elementRef")
                if citation_id and citation_id not in self.citation_ids and citation_id not in self.all_element_ids:
                    self.report.add_error(
                        "MISSING_CITATION_REF",
                        f"Interaction '{interaction_id}' references missing citation '{citation_id}'",
                        interaction_id
                    )

            # Validate Point elementRef references (must reference DataNodes)
            for point in interaction.findall("gpml:Waypoints/gpml:Point", NS):
                point_elem_ref = point.get("elementRef")
                if point_elem_ref and point_elem_ref not in self.datanode_ids and point_elem_ref not in self.all_element_ids:
                    self.report.add_error(
                        "MISSING_DATANODE_REF",
                        f"Point in Interaction '{interaction_id}' references missing DataNode '{point_elem_ref}'",
                        interaction_id
                    )

        # Validate Group groupRef (parent group references)
        for group in root.findall("gpml:Groups/gpml:Group", NS):
            group_id = group.get("elementId", "UNKNOWN")
            parent_group_ref = group.get("groupRef")
            if parent_group_ref and parent_group_ref not in self.group_ids and parent_group_ref not in self.all_element_ids:
                self.report.add_error(
                    "MISSING_PARENT_GROUP_REF",
                    f"Group '{group_id}' references missing parent group '{parent_group_ref}'",
                    group_id
                )

            # Validate CitationRef references
            for citation_ref in group.findall("gpml:CitationRef", NS):
                citation_id = citation_ref.get("elementRef")
                if citation_id and citation_id not in self.citation_ids and citation_id not in self.all_element_ids:
                    self.report.add_error(
                        "MISSING_CITATION_REF",
                        f"Group '{group_id}' references missing citation '{citation_id}'",
                        group_id
                    )

    def _validate_required_fields(self, root: ET.Element):
        """Validate required fields for each element type"""

        # DataNode validation
        for datanode in root.findall("gpml:DataNodes/gpml:DataNode", NS):
            datanode_id = datanode.get("elementId")
            if not datanode_id:
                self.report.add_error("MISSING_REQUIRED_FIELD", "DataNode missing 'elementId'")

            text_label = datanode.get("textLabel")
            if not text_label:
                self.report.add_error(
                    "MISSING_REQUIRED_FIELD",
                    f"DataNode '{datanode_id}' missing 'textLabel'",
                    datanode_id
                )

            node_type = datanode.get("type")
            if not node_type:
                self.report.add_error(
                    "MISSING_REQUIRED_FIELD",
                    f"DataNode '{datanode_id}' missing 'type'",
                    datanode_id
                )
            elif node_type not in VALID_DATANODE_TYPES:
                self.report.add_warning(
                    f"DataNode '{datanode_id}' has unknown type '{node_type}'"
                )

            # Check for graphics
            graphics = datanode.find("gpml:Graphics", NS)
            if graphics is None:
                self.report.add_error(
                    "MISSING_REQUIRED_FIELD",
                    f"DataNode '{datanode_id}' missing Graphics element",
                    datanode_id
                )
            else:
                for attr in ["centerX", "centerY", "width", "height"]:
                    if graphics.get(attr) is None:
                        self.report.add_error(
                            "MISSING_REQUIRED_FIELD",
                            f"DataNode '{datanode_id}' Graphics missing '{attr}'",
                            datanode_id
                        )

        # Interaction validation
        for interaction in root.findall("gpml:Interactions/gpml:Interaction", NS):
            interaction_id = interaction.get("elementId")
            if not interaction_id:
                self.report.add_error("MISSING_REQUIRED_FIELD", "Interaction missing 'elementId'")

            waypoints = interaction.find("gpml:Waypoints", NS)
            if waypoints is None:
                self.report.add_error(
                    "MISSING_REQUIRED_FIELD",
                    f"Interaction '{interaction_id}' missing Waypoints",
                    interaction_id
                )
            else:
                points = waypoints.findall("gpml:Point", NS)
                if len(points) < 2:
                    self.report.add_error(
                        "MISSING_REQUIRED_FIELD",
                        f"Interaction '{interaction_id}' has fewer than 2 waypoints ({len(points)})",
                        interaction_id
                    )

            graphics = interaction.find("gpml:Graphics", NS)
            if graphics is None:
                self.report.add_error(
                    "MISSING_REQUIRED_FIELD",
                    f"Interaction '{interaction_id}' missing Graphics element",
                    interaction_id
                )
            else:
                for attr in ["lineColor", "lineStyle", "lineWidth", "connectorType"]:
                    if graphics.get(attr) is None:
                        self.report.add_error(
                            "MISSING_REQUIRED_FIELD",
                            f"Interaction '{interaction_id}' Graphics missing '{attr}'",
                            interaction_id
                        )

        # Group validation
        for group in root.findall("gpml:Groups/gpml:Group", NS):
            group_id = group.get("elementId")
            if not group_id:
                self.report.add_error("MISSING_REQUIRED_FIELD", "Group missing 'elementId'")

            graphics = group.find("gpml:Graphics", NS)
            if graphics is None:
                self.report.add_error(
                    "MISSING_REQUIRED_FIELD",
                    f"Group '{group_id}' missing Graphics element",
                    group_id
                )
            else:
                for attr in ["centerX", "centerY", "width", "height"]:
                    if graphics.get(attr) is None:
                        self.report.add_error(
                            "MISSING_REQUIRED_FIELD",
                            f"Group '{group_id}' Graphics missing '{attr}'",
                            group_id
                        )

        # Citation validation
        for citation in root.findall("gpml:Citations/gpml:Citation", NS):
            citation_id = citation.get("elementId")
            if not citation_id:
                self.report.add_error("MISSING_REQUIRED_FIELD", "Citation missing 'elementId'")

            xref = citation.find("gpml:Xref", NS)
            url = citation.find("gpml:Url", NS)
            if xref is None and url is None:
                self.report.add_error(
                    "MISSING_REQUIRED_FIELD",
                    f"Citation '{citation_id}' must have either Xref or Url",
                    citation_id
                )

    def _validate_graphics(self, root: ET.Element):
        """Validate graphics attributes for valid values"""

        # Validate Point arrow heads
        for point in root.findall(".//gpml:Point", NS):
            arrow_head = point.get("arrowHead")
            if arrow_head and arrow_head not in VALID_ARROW_HEADS:
                self.report.add_warning(
                    f"Unknown arrowHead value: '{arrow_head}'. Expected one of: {VALID_ARROW_HEADS}"
                )

        # Validate coordinate values are numeric
        for element in root.findall(".//gpml:Graphics", NS):
            for attr in ["centerX", "centerY", "width", "height", "x", "y"]:
                value = element.get(attr)
                if value:
                    try:
                        float(value)
                    except ValueError:
                        parent = element.getparent() if hasattr(element, 'getparent') else None
                        self.report.add_warning(
                            f"Graphics attribute '{attr}' has non-numeric value: '{value}'"
                        )

    def _validate_xrefs(self, root: ET.Element):
        """Validate cross-references"""
        for element in root.findall(".//gpml:Xref", NS):
            parent = element
            # Navigate up to find parent's elementId
            elem_id = "unknown"

            # Check for required attributes
            identifier = element.get("identifier")
            if not identifier:
                self.report.add_warning("Xref missing 'identifier' attribute")

            data_source = element.get("dataSource")
            if not data_source:
                self.report.add_warning("Xref missing 'dataSource' attribute")

    def _validate_circular_references(self, root: ET.Element):
        """Check for circular groupRef references"""
        # Build a graph of group references
        group_graph = {}
        for group in root.findall("gpml:Groups/gpml:Group", NS):
            group_id = group.get("elementId")
            parent_ref = group.get("groupRef")
            if group_id:
                group_graph[group_id] = parent_ref

        # Check for cycles using DFS
        def has_cycle(node, visited, rec_stack):
            visited.add(node)
            rec_stack.add(node)

            # Visit parent
            parent = group_graph.get(node)
            if parent and parent in group_graph:
                if parent not in visited:
                    if has_cycle(parent, visited, rec_stack):
                        return True
                elif parent in rec_stack:
                    return True

            rec_stack.remove(node)
            return False

        visited = set()
        for group_id in group_graph:
            if group_id not in visited:
                if has_cycle(group_id, visited, set()):
                    self.report.add_error(
                        "CIRCULAR_GROUP_REF",
                        f"Circular groupRef reference detected involving group '{group_id}'",
                        group_id
                    )

    def _validate_missing_elements(self, root: ET.Element):
        """Check for elements with missing required elementId attributes"""
        # Check DataNodes
        for i, datanode in enumerate(root.findall("gpml:DataNodes/gpml:DataNode", NS)):
            if not datanode.get("elementId"):
                text_label = datanode.get("textLabel", "UNKNOWN")
                self.report.add_error(
                    "MISSING_ELEMENT_ID",
                    f"DataNode #{i+1} ('{text_label}') missing elementId attribute"
                )

        # Check Interactions
        for i, interaction in enumerate(root.findall("gpml:Interactions/gpml:Interaction", NS)):
            if not interaction.get("elementId"):
                self.report.add_error(
                    "MISSING_ELEMENT_ID",
                    f"Interaction #{i+1} missing elementId attribute"
                )

        # Check Groups
        for i, group in enumerate(root.findall("gpml:Groups/gpml:Group", NS)):
            if not group.get("elementId"):
                text_label = group.get("textLabel", "UNKNOWN")
                self.report.add_error(
                    "MISSING_ELEMENT_ID",
                    f"Group #{i+1} ('{text_label}') missing elementId attribute"
                )

        # Check Citations
        for i, citation in enumerate(root.findall("gpml:Citations/gpml:Citation", NS)):
            if not citation.get("elementId"):
                self.report.add_error(
                    "MISSING_ELEMENT_ID",
                    f"Citation #{i+1} missing elementId attribute"
                )


def validate_directory(directory: str) -> List[ValidationReport]:
    """Validate all GPML files in a directory"""
    reports = []
    gpml_files = list(Path(directory).rglob("*.gpml"))

    if not gpml_files:
        print(f"No GPML files found in {directory}")
        return reports

    validator = GPMLValidator()
    for gpml_file in sorted(gpml_files):
        print(f"Validating: {gpml_file}")
        report = validator.validate_file(str(gpml_file))
        reports.append(report)

    return reports


def print_report(report: ValidationReport):
    """Print a validation report"""
    print(f"\n{'='*80}")
    print(f"File: {report.file_path}")
    print(f"{'='*80}")

    # Statistics
    print(f"\nStatistics:")
    for key, value in report.stats.items():
        print(f"  {key}: {value}")

    # XML validity
    print(f"\nXML Valid: {report.xml_valid}")

    # Errors
    if report.errors:
        print(f"\nErrors ({len(report.errors)}):")
        error_types = defaultdict(int)
        for error in report.errors:
            error_types[error.error_type] += 1

        for error_type, count in sorted(error_types.items()):
            print(f"  {error_type}: {count}")

        # Show detailed errors
        if len(report.errors) <= 200:
            print("\nDetailed Errors:")
            for error in report.errors:
                elem_info = f" (ID: {error.element_id})" if error.element_id else ""
                print(f"  - [{error.error_type}]{elem_info}: {error.message}")
        else:
            print(f"\nShowing first 20 of {len(report.errors)} errors:")
            for error in report.errors[:20]:
                elem_info = f" (ID: {error.element_id})" if error.element_id else ""
                print(f"  - [{error.error_type}]{elem_info}: {error.message}")
    else:
        print("\nNo errors found!")

    # Warnings
    if report.warnings:
        print(f"\nWarnings ({len(report.warnings)}):")
        for warning in report.warnings[:200]:
            print(f"  - {warning}")
        if len(report.warnings) > 200:
            print(f"  ... and {len(report.warnings) - 200} more warnings")

    # Summary
    status = "VALID" if report.is_valid() else "INVALID"
    print(f"\nStatus: {status}\n")


def print_summary(reports: List[ValidationReport]):
    """Print summary of all reports"""
    print(f"\n{'='*80}")
    print("VALIDATION SUMMARY")
    print(f"{'='*80}\n")

    total_files = len(reports)
    valid_files = sum(1 for r in reports if r.is_valid())
    total_errors = sum(len(r.errors) for r in reports)
    total_warnings = sum(len(r.warnings) for r in reports)

    print(f"Total files: {total_files}")
    print(f"Valid files: {valid_files}")
    print(f"Invalid files: {total_files - valid_files}")
    print(f"Total errors: {total_errors}")
    print(f"Total warnings: {total_warnings}")

    # Error breakdown
    if total_errors > 0:
        all_error_types = defaultdict(int)
        for report in reports:
            for error in report.errors:
                all_error_types[error.error_type] += 1

        print(f"\nError Types:")
        for error_type, count in sorted(all_error_types.items(), key=lambda x: -x[1]):
            print(f"  {error_type}: {count}")

    print()


def main():
    """Main entry point"""
    if len(sys.argv) < 2:
        print("Usage: python test_gpml_files.py <gpml_file_or_directory>")
        sys.exit(1)

    path = sys.argv[1]

    if os.path.isfile(path):
        # Validate single file
        validator = GPMLValidator()
        report = validator.validate_file(path)
        print_report(report)
        sys.exit(0 if report.is_valid() else 1)

    elif os.path.isdir(path):
        # Validate directory
        reports = validate_directory(path)

        # Print reports
        for report in reports:
            print_report(report)

        # Print summary
        print_summary(reports)

        # Exit with error if any invalid files
        any_invalid = any(not r.is_valid() for r in reports)
        sys.exit(1 if any_invalid else 0)

    else:
        print(f"Error: '{path}' is not a file or directory")
        sys.exit(1)


if __name__ == "__main__":
    main()
