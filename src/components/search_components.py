"""
Shared search component utilities for consistent autocomplete behavior across pages.

This module provides reusable autocomplete components and helper functions
for creating consistent search interfaces in wiki, synteny, and other pages.
"""

import dash_mantine_components as dmc
from dash_iconify import DashIconify
import pandas as pd
from typing import List, Dict, Optional


def create_search_autocomplete(
    id: str,
    label: str,
    placeholder: str,
    limit: int = 20,
    icon: str = "bi:search",
    **kwargs,
) -> dmc.Autocomplete:
    """
    Create a standardized autocomplete search component.

    Args:
        id: Component ID for callbacks
        label: Label text displayed above the input
        placeholder: Placeholder text shown when empty
        limit: Maximum number of dropdown items to show (default: 20)
        icon: Icon to display on the left (default: "bi:search")
        **kwargs: Additional props to pass to dmc.Autocomplete

    Returns:
        dmc.Autocomplete component configured with consistent styling

    Example:
        >>> taxonomy_search = create_search_autocomplete(
        ...     id="taxonomy-autocomplete",
        ...     label="Search Taxonomy",
        ...     placeholder="Type to search species, genus, family...",
        ... )
    """
    default_props = {
        "data": [],  # Will be populated by callbacks
        "limit": limit,
        "style": {"width": "100%"},
        "leftSection": DashIconify(icon=icon),
    }

    # Merge default props with user-provided kwargs
    props = {**default_props, **kwargs}

    return dmc.Autocomplete(id=id, label=label, placeholder=placeholder, **props)


def get_unique_values_from_columns(
    df: pd.DataFrame,
    columns: List[str],
    sort_by_length: bool = True,
    max_length: Optional[int] = None,
) -> List[Dict[str, str]]:
    """
    Extract unique non-null values from specified DataFrame columns.

    Args:
        df: DataFrame to extract values from
        columns: List of column names to search
        sort_by_length: If True, sort by string length (shorter first)
        max_length: Optional maximum string length to include

    Returns:
        List of dicts with 'value' and 'label' keys for autocomplete data

    Example:
        >>> df = pd.DataFrame({'genus': ['Fusarium', 'Aspergillus'], 'family': ['Nectriaceae', 'Aspergillaceae']})
        >>> get_unique_values_from_columns(df, ['genus', 'family'])
        [{'value': 'Fusarium', 'label': 'Fusarium'}, ...]
    """
    all_values = set()

    for col in columns:
        if col in df.columns:
            # Get non-null values and add to set
            values = df[col].dropna().astype(str).unique()

            # Filter by length if specified
            if max_length:
                values = [v for v in values if len(v) <= max_length]

            all_values.update(values)

    # Convert to sorted list
    if sort_by_length:
        sorted_values = sorted(all_values, key=lambda s: len(s))
    else:
        sorted_values = sorted(all_values)

    # Format for Autocomplete component
    return [{"value": val, "label": val} for val in sorted_values]


def get_unique_values_from_dict_list(
    data: List[Dict],
    keys: List[str],
    sort_by_length: bool = True,
    max_length: Optional[int] = None,
    exclude_empty: bool = True,
) -> List[Dict[str, str]]:
    """
    Extract unique non-empty values from specified keys in a list of dictionaries.

    Useful when working with data directly from database queries (not DataFrames).

    Args:
        data: List of dictionaries
        keys: List of dictionary keys to extract values from
        sort_by_length: If True, sort by string length (shorter first)
        max_length: Optional maximum string length to include
        exclude_empty: If True, exclude empty strings

    Returns:
        List of dicts with 'value' and 'label' keys for autocomplete data

    Example:
        >>> data = [{'name': 'Fusarium', 'family': 'Voyager'}, {'name': 'Aspergillus', 'family': 'Atlantia'}]
        >>> get_unique_values_from_dict_list(data, ['name', 'family'])
        [{'value': 'Voyager', 'label': 'Voyager'}, ...]
    """
    all_values = set()

    for item in data:
        for key in keys:
            value = item.get(key)

            # Skip None and optionally empty strings
            if value is None:
                continue
            if exclude_empty and str(value).strip() == "":
                continue

            str_value = str(value)

            # Filter by length if specified
            if max_length and len(str_value) > max_length:
                continue

            all_values.add(str_value)

    # Convert to sorted list
    if sort_by_length:
        sorted_values = sorted(all_values, key=lambda s: len(s))
    else:
        sorted_values = sorted(all_values)

    # Format for Autocomplete component
    return [{"value": val, "label": val} for val in sorted_values]


def filter_dataframe_by_search(
    df: pd.DataFrame,
    search_value: Optional[str],
    columns: List[str],
    case_sensitive: bool = False,
    exact_match: bool = False,
) -> pd.DataFrame:
    """
    Filter a DataFrame based on search value across multiple columns.

    Args:
        df: DataFrame to filter
        search_value: Value to search for (None or empty string returns unfiltered df)
        columns: List of column names to search in
        case_sensitive: If True, perform case-sensitive matching
        exact_match: If True, require exact match; otherwise partial match

    Returns:
        Filtered DataFrame

    Example:
        >>> df = pd.DataFrame({'name': ['Fusarium graminearum', 'Aspergillus niger']})
        >>> filtered = filter_dataframe_by_search(df, 'Fusarium', ['name'])
        >>> len(filtered)
        1
    """
    # Return original if no search value
    if not search_value or not search_value.strip():
        return df

    # Create mask for rows matching search value
    mask = pd.Series([False] * len(df))

    for col in columns:
        if col not in df.columns:
            continue

        if exact_match:
            # Exact matching
            if case_sensitive:
                mask |= df[col].astype(str) == search_value
            else:
                mask |= df[col].astype(str).str.lower() == search_value.lower()
        else:
            # Partial matching (contains)
            mask |= (
                df[col]
                .astype(str)
                .str.contains(
                    search_value,
                    case=case_sensitive,
                    na=False,
                    regex=False,  # Treat search value as literal string, not regex
                )
            )

    return df[mask]


def filter_dict_list_by_search(
    data: List[Dict],
    search_value: Optional[str],
    keys: List[str],
    case_sensitive: bool = False,
    exact_match: bool = False,
) -> List[Dict]:
    """
    Filter a list of dictionaries based on search value across multiple keys.

    Args:
        data: List of dictionaries to filter
        search_value: Value to search for (None or empty string returns original list)
        keys: List of dictionary keys to search in
        case_sensitive: If True, perform case-sensitive matching
        exact_match: If True, require exact match; otherwise partial match

    Returns:
        Filtered list of dictionaries

    Example:
        >>> data = [{'name': 'Fusarium', 'family': 'Voyager'}, {'name': 'Aspergillus', 'family': 'Atlantia'}]
        >>> filtered = filter_dict_list_by_search(data, 'Voyager', ['family'])
        >>> len(filtered)
        1
    """
    # Return original if no search value
    if not search_value or not search_value.strip():
        return data

    search_term = search_value if case_sensitive else search_value.lower()
    filtered = []

    for item in data:
        match_found = False

        for key in keys:
            value = item.get(key)
            if value is None:
                continue

            str_value = str(value) if case_sensitive else str(value).lower()

            if exact_match:
                if str_value == search_term:
                    match_found = True
                    break
            else:
                if search_term in str_value:
                    match_found = True
                    break

        if match_found:
            filtered.append(item)

    return filtered
