"""
Text Label Cleaner Module

This module provides functions to clean and format text labels for GPML output,
converting HTML-like tags and entities to proper Unicode characters.
"""

import re
import html

class TextLabelCleaner:
    """
    Cleans and formats text labels for better display in GPML.
    Converts HTML entities, superscripts, subscripts, and other formatting.
    """

    # Unicode superscript mapping
    SUPERSCRIPT_MAP = {
        '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴',
        '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹',
        '+': '⁺', '-': '⁻', '=': '⁼', '(': '⁽', ')': '⁾',
        'n': 'ⁿ', 'i': 'ⁱ'
    }

    # Unicode subscript mapping
    SUBSCRIPT_MAP = {
        '0': '₀', '1': '₁', '2': '₂', '3': '₃', '4': '₄',
        '5': '₅', '6': '₆', '7': '₇', '8': '₈', '9': '₉',
        '+': '₊', '-': '₋', '=': '₌', '(': '₍', ')': '₎',
        'a': 'ₐ', 'e': 'ₑ', 'o': 'ₒ', 'x': 'ₓ', 'h': 'ₕ',
        'k': 'ₖ', 'l': 'ₗ', 'm': 'ₘ', 'n': 'ₙ', 'p': 'ₚ',
        's': 'ₛ', 't': 'ₜ'
    }

    # Greek letter mapping
    GREEK_LETTERS = {
        'alpha': 'α', 'Alpha': 'Α',
        'beta': 'β', 'Beta': 'Β',
        'gamma': 'γ', 'Gamma': 'Γ',
        'delta': 'δ', 'Delta': 'Δ',
        'epsilon': 'ε', 'Epsilon': 'Ε',
        'zeta': 'ζ', 'Zeta': 'Ζ',
        'eta': 'η', 'Eta': 'Η',
        'theta': 'θ', 'Theta': 'Θ',
        'iota': 'ι', 'Iota': 'Ι',
        'kappa': 'κ', 'Kappa': 'Κ',
        'lambda': 'λ', 'Lambda': 'Λ',
        'mu': 'μ', 'Mu': 'Μ',
        'nu': 'ν', 'Nu': 'Ν',
        'xi': 'ξ', 'Xi': 'Ξ',
        'omicron': 'ο', 'Omicron': 'Ο',
        'pi': 'π', 'Pi': 'Π',
        'rho': 'ρ', 'Rho': 'Ρ',
        'sigma': 'σ', 'Sigma': 'Σ',
        'tau': 'τ', 'Tau': 'Τ',
        'upsilon': 'υ', 'Upsilon': 'Υ',
        'phi': 'φ', 'Phi': 'Φ',
        'chi': 'χ', 'Chi': 'Χ',
        'psi': 'ψ', 'Psi': 'Ψ',
        'omega': 'ω', 'Omega': 'Ω'
    }

    # Common HTML entities in biochemistry
    SPECIAL_ENTITIES = {
        '&alpha;': 'α', '&Alpha;': 'Α',
        '&beta;': 'β', '&Beta;': 'Β',
        '&gamma;': 'γ', '&Gamma;': 'Γ',
        '&delta;': 'δ', '&Delta;': 'Δ',
        '&epsilon;': 'ε', '&Epsilon;': 'Ε',
        '&omega;': 'ω', '&Omega;': 'Ω',
        '&mu;': 'μ', '&Mu;': 'Μ',
        '&pi;': 'π', '&Pi;': 'Π',
        '&sigma;': 'σ', '&Sigma;': 'Σ',
        '&tau;': 'τ', '&Tau;': 'Τ',
        '&rarr;': '→', '&larr;': '←', '&harr;': '↔',
        '&uarr;': '↑', '&darr;': '↓',
        '&plusmn;': '±', '&times;': '×',
        '&deg;': '°', '&prime;': '′', '&Prime;': '″',
        '&infin;': '∞', '&asymp;': '≈',
        '&ne;': '≠', '&le;': '≤', '&ge;': '≥',
        '&middot;': '·', '&bull;': '•'
    }

    @classmethod
    def clean_label(cls, text):
        """
        Clean and format a text label for GPML display.

        Args:
            text (str): Raw text label with possible HTML tags and entities

        Returns:
            str: Cleaned text with Unicode formatting
        """
        if not text:
            return text

        # Convert to string if not already
        text = str(text)

        # Step 1: Convert HTML entities
        text = cls._convert_html_entities(text)

        # Step 2: Convert superscript tags
        text = cls._convert_superscripts(text)

        # Step 3: Convert subscript tags
        text = cls._convert_subscripts(text)

        # Step 4: Handle italic tags (remove or replace with markers)
        text = cls._handle_italic_tags(text)

        # Step 5: Handle bold tags
        text = cls._handle_bold_tags(text)

        # Step 6: Clean up any remaining HTML tags
        text = cls._remove_remaining_html_tags(text)

        # Step 7: Fix common chemical notation
        text = cls._fix_chemical_notation(text)

        # Step 8: Clean up whitespace
        text = ' '.join(text.split())

        return text

    @classmethod
    def _convert_html_entities(cls, text):
        """Convert HTML entities to Unicode characters."""
        # First convert known special entities
        for entity, unicode_char in cls.SPECIAL_ENTITIES.items():
            text = text.replace(entity, unicode_char)

        # Then use html.unescape for any remaining standard entities
        text = html.unescape(text)

        return text

    @classmethod
    def _convert_superscripts(cls, text):
        """Convert <SUP> and <sup> tags to Unicode superscripts."""
        # Find all superscript tags
        pattern = re.compile(r'<(?:SUP|sup)>(.*?)</(?:SUP|sup)>', re.IGNORECASE)

        def replace_superscript(match):
            content = match.group(1)
            result = ''
            for char in content:
                if char in cls.SUPERSCRIPT_MAP:
                    result += cls.SUPERSCRIPT_MAP[char]
                else:
                    # If no superscript available, use regular character
                    result += char
            return result

        text = pattern.sub(replace_superscript, text)
        return text

    @classmethod
    def _convert_subscripts(cls, text):
        """Convert <SUB> and <sub> tags to Unicode subscripts."""
        # Find all subscript tags
        pattern = re.compile(r'<(?:SUB|sub)>(.*?)</(?:SUB|sub)>', re.IGNORECASE)

        def replace_subscript(match):
            content = match.group(1)
            result = ''
            for char in content:
                if char in cls.SUBSCRIPT_MAP:
                    result += cls.SUBSCRIPT_MAP[char]
                else:
                    # If no subscript available, use regular character
                    result += char
            return result

        text = pattern.sub(replace_subscript, text)
        return text

    @classmethod
    def _handle_italic_tags(cls, text):
        """Handle italic tags - remove them but keep the content."""
        # Remove <i> and <em> tags but keep content
        text = re.sub(r'<(?:i|I|em|EM)>', '', text)
        text = re.sub(r'</(?:i|I|em|EM)>', '', text)
        return text

    @classmethod
    def _handle_bold_tags(cls, text):
        """Handle bold tags - remove them but keep the content."""
        # Remove <b> and <strong> tags but keep content
        text = re.sub(r'<(?:b|B|strong|STRONG)>', '', text)
        text = re.sub(r'</(?:b|B|strong|STRONG)>', '', text)
        return text

    @classmethod
    def _remove_remaining_html_tags(cls, text):
        """Remove any remaining HTML tags."""
        # Remove any remaining tags
        text = re.sub(r'<[^>]+>', '', text)
        return text

    @classmethod
    def _fix_chemical_notation(cls, text):
        """Fix common chemical notation patterns."""
        # Convert common patterns like (R)- or (S)- to use proper formatting
        text = re.sub(r'\(R\)-', '(R)-', text)
        text = re.sub(r'\(S\)-', '(S)-', text)
        text = re.sub(r'\(L\)-', 'L-', text)
        text = re.sub(r'\(D\)-', 'D-', text)

        # Fix spacing around operators
        text = re.sub(r'\s*\+\s*', '+', text)  # Remove spaces around +
        text = re.sub(r'\s*-\s*', '-', text)  # Remove spaces around -

        # Fix common chemical group notation
        text = text.replace('orthophosphate', 'phosphate')
        text = text.replace('diphosphate', 'pyrophosphate')

        return text

    @classmethod
    def clean_reaction_equation(cls, equation):
        """
        Clean a reaction equation for display.

        Args:
            equation (str): Raw reaction equation

        Returns:
            str: Cleaned equation with proper arrows
        """
        if not equation:
            return equation

        # Convert arrow notations
        equation = equation.replace('<=>', '⇌')  # Reversible
        equation = equation.replace('<->', '⇌')  # Reversible
        equation = equation.replace('->', '→')  # Forward
        equation = equation.replace('=>', '→')  # Forward
        equation = equation.replace('<-', '←')  # Reverse
        equation = equation.replace('<=', '←')  # Reverse

        # Clean compound names in the equation
        parts = re.split(r'([⇌→←])', equation)
        cleaned_parts = []

        for part in parts:
            if part in ['⇌', '→', '←']:
                cleaned_parts.append(part)
            else:
                # This is a compound list
                compounds = part.split('+')
                cleaned_compounds = []
                for compound in compounds:
                    compound = compound.strip()
                    if compound:
                        # Check for coefficient
                        match = re.match(r'^(\d+(?:\.\d+)?)\s*(.+)$', compound)
                        if match:
                            coeff = match.group(1)
                            name = cls.clean_label(match.group(2))
                            cleaned_compounds.append(f"{coeff} {name}")
                        else:
                            cleaned_compounds.append(cls.clean_label(compound))
                cleaned_parts.append(' + '.join(cleaned_compounds))

        return ' '.join(cleaned_parts)


# Convenience function for easy import
def clean_text_label(text):
    """
    Clean a text label for GPML display.

    Args:
        text: Raw text label

    Returns:
        str: Cleaned text
    """
    return TextLabelCleaner.clean_label(text)


def clean_reaction_equation(equation):
    """
    Clean a reaction equation for display.

    Args:
        equation: Raw reaction equation

    Returns:
        str: Cleaned equation
    """
    return TextLabelCleaner.clean_reaction_equation(equation)
