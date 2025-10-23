import chardet
from pathlib import Path
from typing import Optional, List, Dict, Any


class FileReader:
    """Find and read a file with encoding handling."""

    def __init__(self, base_dir: Optional[Path] = None):
        if base_dir is None:
            base_dir = Path(__file__).resolve().parent.parent.parent
        self.base_dir = Path(base_dir)

    def find(self, filename: str) -> Optional[Path]:
        """Search for the file recursively under base_dir."""
        matches = list(self.base_dir.rglob(filename))
        if not matches:
            print(f"File '{filename}' not found under {self.base_dir}")
            return None
        print(f"Found file: {matches[0]}")
        return matches[0]

    def _detect_encoding(self, file_path: Path) -> Optional[str]:
        with open(file_path, 'rb') as f:
            raw_data = f.read(10000)
            result = chardet.detect(raw_data)
        encoding = result['encoding']
        confidence = result['confidence']
        print(f"Detected encoding for {file_path.name}: {encoding} (confidence: {confidence:.2f})")
        return encoding if confidence >= 0.7 else None

    def read(self, filename: str) -> str:
        """Find the file and read it with the appropriate encoding."""
        file_path = self.find(filename)
        if file_path is None:
            raise FileNotFoundError(f"Could not find file: {filename}")

        detected_encoding = self._detect_encoding(file_path)
        # For enzrxns.dat, prioritize latin-1 encoding
        if 'enzrxns' in filename.lower():
            encodings_to_try = ['latin-1', 'iso-8859-1', 'cp1252', 'utf-8']
        else:
            encodings_to_try = [e for e in [detected_encoding, 'latin-1', 'utf-8', 'iso-8859-1', 'cp1252'] if e]
            encodings_to_try = list(dict.fromkeys(encodings_to_try))

        for enc in encodings_to_try:
            try:
                with open(file_path, 'r', encoding=enc) as f:
                    content = f.read()
                print(f"Successfully read {file_path.name} with encoding: {enc}")
                return content
            except Exception as e:
                if 'enzrxns' not in filename.lower():
                    print(f"Failed with {enc}: {e}")
                continue

        print(f"All encodings failed, trying utf-8 with error replacement...")
        with open(file_path, 'r', encoding='utf-8', errors='replace') as f:
            return f.read()


class RecordParser:
    """Parse structured records from file content."""

    def __init__(self, record_separator: str = "\n//", comment_prefix: str = "#"):
        self.record_separator = record_separator
        self.comment_prefix = comment_prefix

    def split_records(self, content: str) -> List[str]:
        """Split content into individual records, excluding header comments."""
        parts = content.split(self.record_separator)

        records = []
        for part in parts:
            part = part.strip()
            if not part:
                continue

            lines = part.split('\n')
            non_comment_lines = [line for line in lines if not line.strip().startswith(self.comment_prefix)]

            if non_comment_lines:
                clean_record = '\n'.join(non_comment_lines).strip()
                if clean_record:
                    records.append(clean_record)

        return records

    def _fix_split_words(self, text: str) -> str:
        """
        Fix words that are incorrectly split across lines with newlines.

        Simple approach: Replace all newlines with spaces, then normalize whitespace.
        This preserves word boundaries correctly.

        Examples:
            "T\nhe" -> "T he" (will be caught by normalization)
            "terpenoid\nindole" -> "terpenoid indole"
            "families\nApocynaceae" -> "families Apocynaceae"
            "terpenoid \nindole" -> "terpenoid indole"
        """
        import re

        # Step 1: Replace all newlines with spaces
        text = text.replace('\n', ' ')

        # Step 2: Collapse multiple spaces into one
        text = re.sub(r'\s+', ' ', text)

        # Step 3: Strip leading/trailing whitespace
        text = text.strip()

        return text

    def parse_record(self, record_text: str) -> Dict[str, Any]:
        """Parse a single record into a dictionary, handling multi-line values with / continuation."""
        record = {}
        lines = record_text.strip().split('\n')

        current_key = None
        current_value = []

        for line in lines:
            original_line = line
            line = line.strip()
            if not line or line.startswith(self.comment_prefix):
                continue

            # Check if this is a new field line (has ' - ' and doesn't start with /)
            is_new_field = ' - ' in line and not line.startswith('/')

            # If we encounter a new field line, save the previous key-value
            if is_new_field and current_key:
                value = '\n'.join(current_value)
                # Fix split words in text fields
                if current_key in ['COMMENT', 'DESCRIPTION', 'CREDITS']:
                    value = self._fix_split_words(value)
                if current_key in record:
                    if not isinstance(record[current_key], list):
                        record[current_key] = [record[current_key]]
                    record[current_key].append(value)
                else:
                    record[current_key] = value
                current_key = None
                current_value = []

            # Parse new field line
            if is_new_field:
                key, value = line.split(' - ', 1)
                current_key = key.strip()
                current_value = [value.strip()]
            # Check if this is a continuation line (starts with / and we have a current key)
            elif line.startswith('/') and current_key:
                # Append continuation to current value (remove leading /)
                current_value.append(line[1:])

        # Don't forget the last key-value pair
        if current_key:
            value = '\n'.join(current_value)
            # Fix split words in text fields
            if current_key in ['COMMENT', 'DESCRIPTION', 'CREDITS']:
                value = self._fix_split_words(value)
            if current_key in record:
                if not isinstance(record[current_key], list):
                    record[current_key] = [record[current_key]]
                record[current_key].append(value)
            else:
                record[current_key] = value

        return record

    def parse_all_records(self, content: str) -> List[Dict[str, Any]]:
        """Parse all records from content into a list of dictionaries."""
        records = self.split_records(content)
        return [self.parse_record(record) for record in records]


class DataProcessor:
    """Process and analyze parsed records."""

    def __init__(self, records: List[Dict[str, Any]]):
        self.records = records

    def get_unique_keys(self) -> set:
        """Get all unique keys across all records."""
        keys = set()
        for record in self.records:
            keys.update(record.keys())
        return keys

    def filter_by_key(self, key: str, value: str = None) -> List[Dict[str, Any]]:
        """Filter records that have a specific key, optionally with a specific value."""
        filtered = []
        for record in self.records:
            if key in record:
                if value is None:
                    filtered.append(record)
                else:
                    record_values = record[key] if isinstance(record[key], list) else [record[key]]
                    if value in record_values:
                        filtered.append(record)
        return filtered

    def get_key_statistics(self) -> Dict[str, int]:
        """Get statistics on how often each key appears."""
        key_counts = {}
        for record in self.records:
            for key in record.keys():
                key_counts[key] = key_counts.get(key, 0) + 1
        return dict(sorted(key_counts.items(), key=lambda x: x[1], reverse=True))

    def find_by_unique_id(self, unique_id: str) -> Optional[Dict[str, Any]]:
        """Find a record by its UNIQUE-ID."""
        for record in self.records:
            if record.get('UNIQUE-ID') == unique_id:
                return record
        return None

    def summary(self) -> Dict[str, Any]:
        """Get a summary of the dataset."""
        return {
            'total_records': len(self.records),
            'unique_keys': len(self.get_unique_keys()),
            'key_statistics': self.get_key_statistics(),
            'sample_keys': list(self.get_unique_keys())[:10]
        }


def read_and_parse(filename: str, base_dir: Optional[Path] = None) -> DataProcessor:
    """Read a file and parse all records into a DataProcessor."""
    reader = FileReader(base_dir)
    content = reader.read(filename)
    parser = RecordParser()
    records = parser.parse_all_records(content)
    return DataProcessor(records)