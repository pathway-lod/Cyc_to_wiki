"""
ID Management utilities for WikiPathways GPML generation.
"""
import re

def sanitize_element_id(element_id):
    """
    Sanitize element IDs for GPML compatibility.

    Args:
        element_id (str): Original element ID

    Returns:
        str: Sanitized element ID safe for GPML
    """
    if not element_id:
        return element_id
    
    # replace invalid characters with underscores
    # restrict to alphanumeric and underscores for maximum compatibility with xs:ID
    sanitized = re.sub(r'[^a-zA-Z0-9_]', '_', element_id)
    
    # ensure it starts with a letter or underscore
    if sanitized and not (sanitized[0].isalpha() or sanitized[0] == '_'):
        sanitized = f"_{sanitized}"
        
    return sanitized


class IDManager:
    """
    Manages ID mapping and sanitization for pathway elements.

    Ensures all element IDs are unique and GPML-compatible while maintaining
    a mapping between original and sanitized IDs.
    """

    def __init__(self):
        """Initialize empty ID mapping."""
        self.id_mapping = {}
        self.sanitized_ids = {}  # Track which sanitized IDs are in use: {sanitized_id: original_id}
        self.id_counter = {}  # Counter for handling collisions: {base_id: count}

    def register_id(self, original_id):
        """
        Register and sanitize an ID, ensuring uniqueness.

        Args:
            original_id (str): Original ID from BioCyc data

        Returns:
            str: Sanitized ID safe for GPML (guaranteed unique)
        """
        if original_id not in self.id_mapping:
            sanitized = sanitize_element_id(original_id)

            # Check for collision with existing sanitized IDs
            if sanitized in self.sanitized_ids:
                # Collision detected - append counter to make unique
                base_id = sanitized
                if base_id not in self.id_counter:
                    self.id_counter[base_id] = 1
                else:
                    self.id_counter[base_id] += 1

                # Append counter to make unique
                sanitized = f"{base_id}_{self.id_counter[base_id]}"

            self.id_mapping[original_id] = sanitized
            self.sanitized_ids[sanitized] = original_id

        return self.id_mapping[original_id]

    def get_sanitized_id(self, original_id):
        """
        Get sanitized ID for an original ID.

        Args:
            original_id (str): Original ID to look up

        Returns:
            str: Sanitized ID if found, otherwise original ID
        """
        return self.id_mapping.get(original_id, original_id)

    def get_all_sanitized_ids(self):
        """
        Get all currently registered sanitized IDs.

        Returns:
            set: Set of all sanitized IDs in use
        """
        return set(self.sanitized_ids.keys())