"""
Database query helpers for synteny visualization
This module provides functions to query GFF data with all necessary annotations
"""

from sqlalchemy import func, and_, or_
from src.database.models.schema import (
    Gff,
    Accessions,
    Ships,
    ShipAccessions,
    JoinedShips,
    Taxonomy,
    FamilyNames,
    Navis,
    Haplotype,
)
from src.database.sql_engine import get_starbase_session


def get_all_ships_with_gff():
    """
    Get all ships that have GFF annotation data.
    Returns a list of dictionaries with ship information.
    """
    with get_starbase_session() as session:
        query = (
            session.query(
                Ships.id.label('id'),
                ShipAccessions.ship_accession_display,
                func.max(Taxonomy.name).label('name'),
                func.max(FamilyNames.familyName).label('familyName'),
                func.count(Gff.id).label('gff_count'),
            )
            .join(Gff, Gff.ship_id == Ships.id)
            .outerjoin(ShipAccessions, ShipAccessions.ship_id == Ships.id)
            .outerjoin(JoinedShips, JoinedShips.ship_id == Ships.id)
            .outerjoin(Taxonomy, JoinedShips.tax_id == Taxonomy.id)
            .outerjoin(FamilyNames, JoinedShips.ship_family_id == FamilyNames.id)
            .group_by(Ships.id, ShipAccessions.ship_accession_display)
            .having(func.count(Gff.id) > 0)
            .order_by(func.max(Taxonomy.name), ShipAccessions.ship_accession_display)
        )
        
        results = []
        for row in query.all():
            results.append({
                'id': row.id,
                'ship_accession_display': row.ship_accession_display or '',
                'name': row.name or '',
                'familyName': row.familyName or '',
                'gff_count': row.gff_count,
            })
        
        return results


def get_ships_unified_search():
    """
    Get all ships with GFF data including comprehensive classification information
    for unified search across accessions, taxonomy, family names, navis, and haplotypes.
    
    Returns:
        List of dictionaries with comprehensive ship metadata for search/display
    """
    with get_starbase_session() as session:
        # Subquery to get one joined_ship per ship (to avoid duplicates)
        js_subq = (
            session.query(
                JoinedShips.ship_id.label('ship_id'),
                func.min(JoinedShips.id).label('js_id'),
            )
            .group_by(JoinedShips.ship_id)
            .subquery()
        )
        
        query = (
            session.query(
                Ships.id.label('id'),
                ShipAccessions.ship_accession_display,
                Accessions.accession_tag,
                Accessions.version_tag,
                Taxonomy.name.label('taxonomy_name'),
                FamilyNames.familyName,
                Navis.navis_name,
                Haplotype.haplotype_name,
                func.count(Gff.id).label('gff_count'),
            )
            .join(Gff, Gff.ship_id == Ships.id)
            .outerjoin(ShipAccessions, ShipAccessions.ship_id == Ships.id)
            .outerjoin(Accessions, Ships.accession_id == Accessions.id)
            .outerjoin(js_subq, js_subq.c.ship_id == Ships.id)
            .outerjoin(
                JoinedShips,
                and_(
                    JoinedShips.ship_id == js_subq.c.ship_id,
                    JoinedShips.id == js_subq.c.js_id,
                ),
            )
            .outerjoin(Taxonomy, JoinedShips.tax_id == Taxonomy.id)
            .outerjoin(FamilyNames, JoinedShips.ship_family_id == FamilyNames.id)
            .outerjoin(Navis, JoinedShips.ship_navis_id == Navis.id)
            .outerjoin(Haplotype, JoinedShips.ship_haplotype_id == Haplotype.id)
            .group_by(
                Ships.id,
                ShipAccessions.ship_accession_display,
                Accessions.accession_tag,
                Accessions.version_tag,
                Taxonomy.name,
                FamilyNames.familyName,
                Navis.navis_name,
                Haplotype.haplotype_name,
            )
            .having(func.count(Gff.id) > 0)
            .order_by(Taxonomy.name, FamilyNames.familyName, ShipAccessions.ship_accession_display)
        )
        
        results = []
        for row in query.all():
            # Build accession display
            accession_display = row.ship_accession_display
            if not accession_display and row.accession_tag:
                accession_display = f"{row.accession_tag}.{row.version_tag}"
            
            # Build searchable label: SSB accession first, then taxonomy/family context
            label_parts = []
            if accession_display:
                label_parts.append(accession_display)
            if row.taxonomy_name:
                label_parts.append(f"| {row.taxonomy_name}")
            if row.familyName:
                label_parts.append(f"[{row.familyName}]")
            if row.navis_name:
                label_parts.append(f"({row.navis_name})")
            if row.haplotype_name:
                label_parts.append(f"<{row.haplotype_name}>")

            label = " ".join(label_parts) if label_parts else f"Ship {row.id}"
            
            # Determine grouping for dropdown organization
            group = row.taxonomy_name or "Unknown Taxonomy"
            
            results.append({
                'id': row.id,
                'ship_accession_display': accession_display or '',
                'accession_display': f"{row.accession_tag}.{row.version_tag}" if row.accession_tag else '',
                'taxonomy_name': row.taxonomy_name or '',
                'familyName': row.familyName or '',
                'navis_name': row.navis_name or '',
                'haplotype_name': row.haplotype_name or '',
                'gff_count': row.gff_count,
                'label': label,
                'group': group,
            })
        
        return results


def get_gff_by_ship_ids(ship_ids):
    """
    Get all GFF entries for a list of ship IDs with related accession data.
    
    Args:
        ship_ids: List of ship IDs
    
    Returns:
        List of dictionaries containing GFF data with metadata
    """
    with get_starbase_session() as session:
        # One joined_ship per ship (min id) to avoid duplicating GFF rows
        js_subq = (
            session.query(
                JoinedShips.ship_id.label('ship_id'),
                func.min(JoinedShips.id).label('js_id'),
            )
            .group_by(JoinedShips.ship_id)
            .subquery()
        )
        query = (
            session.query(Gff, Accessions, Ships, ShipAccessions, Taxonomy, FamilyNames)
            .join(Accessions, Gff.accession_id == Accessions.id)
            .join(Ships, Gff.ship_id == Ships.id)
            .outerjoin(ShipAccessions, ShipAccessions.ship_id == Ships.id)
            .outerjoin(
                js_subq,
                js_subq.c.ship_id == Ships.id,
            )
            .outerjoin(
                JoinedShips,
                and_(
                    JoinedShips.ship_id == js_subq.c.ship_id,
                    JoinedShips.id == js_subq.c.js_id,
                ),
            )
            .outerjoin(Taxonomy, JoinedShips.tax_id == Taxonomy.id)
            .outerjoin(FamilyNames, JoinedShips.ship_family_id == FamilyNames.id)
            .filter(Gff.ship_id.in_(ship_ids))
            .order_by(Ships.id, Gff.start)
        )
        
        results = []
        for row in query.all():
            gff, accession, ship, ship_acc, taxonomy, family = row
            results.append({
                # GFF fields
                'id': gff.id,
                'accession_id': gff.accession_id,
                'source': gff.source,
                'type': gff.type,
                'start': gff.start,
                'end': gff.end,
                'phase': gff.phase,
                'strand': gff.strand,
                'score': gff.score,
                'attributes': gff.attributes,
                'ship_id': gff.ship_id,
                # Related accession data (Accessions has accession_tag, not accession)
                'accession': getattr(accession, 'accession_display', None) or (
                    f"{accession.accession_tag}.{accession.version_tag}" if accession else ''
                ),
                'name': taxonomy.name if taxonomy else '',
                # Related ship data from ShipAccessions and FamilyNames
                'ship_accession_display': ship_acc.ship_accession_display if ship_acc else '',
                'familyName': family.familyName if family else '',
            })
        
        return results


def get_gff_by_accession_id(accession_id):
    """
    Get all GFF entries for a specific accession.
    
    Args:
        accession_id: Accession ID
    
    Returns:
        List of dictionaries containing GFF data
    """
    with get_starbase_session() as session:
        query = (
            session.query(Gff)
            .filter(Gff.accession_id == accession_id)
            .order_by(Gff.start)
        )
        
        results = []
        for gff in query.all():
            results.append({
                'id': gff.id,
                'accession_id': gff.accession_id,
                'source': gff.source,
                'type': gff.type,
                'start': gff.start,
                'end': gff.end,
                'phase': gff.phase,
                'strand': gff.strand,
                'score': gff.score,
                'attributes': gff.attributes,
                'ship_id': gff.ship_id,
            })
        
        return results


def parse_gff_attributes(attr_string):
    """
    Parse GFF attributes string into a dictionary.
    
    Args:
        attr_string: Semicolon-separated key=value pairs
    
    Returns:
        Dictionary of attributes
    """
    if not attr_string:
        return {}
    
    attrs = {}
    for item in attr_string.split(';'):
        item = item.strip()
        if '=' in item:
            key, value = item.split('=', 1)
            attrs[key.strip()] = value.strip()
    
    return attrs


def get_gene_families_by_ships(ship_ids):
    """
    Get all unique gene families (Target_IDs) across selected ships.
    Useful for creating color scales.
    
    Args:
        ship_ids: List of ship IDs
    
    Returns:
        List of unique Target_IDs
    """
    gff_data = get_gff_by_ship_ids(ship_ids)
    
    families = set()
    for gff_entry in gff_data:
        attrs = parse_gff_attributes(gff_entry['attributes'])
        target_id = attrs.get('Target_ID')
        if target_id:
            families.add(target_id)
    
    return sorted(list(families))


def get_sequence_statistics(ship_id):
    """
    Get statistics about a genomic sequence.
    
    Args:
        ship_id: Ship ID
    
    Returns:
        Dictionary with statistics
    """
    with get_starbase_session() as session:
        gff_entries = (
            session.query(Gff)
            .filter(Gff.ship_id == ship_id)
            .all()
        )
        
        if not gff_entries:
            return None
        
        starts = [g.start for g in gff_entries]
        ends = [g.end for g in gff_entries]
        types = [g.type for g in gff_entries]
        
        stats = {
            'ship_id': ship_id,
            'num_features': len(gff_entries),
            'min_position': min(starts),
            'max_position': max(ends),
            'total_length': max(ends) - min(starts),
            'feature_types': {t: types.count(t) for t in set(types)},
        }
        
        return stats


def search_genes_by_target(target_id, ship_ids=None):
    """
    Search for genes by Target_ID across ships.
    Useful for finding homologous genes.
    
    Args:
        target_id: Target_ID to search for
        ship_ids: Optional list of ship IDs to limit search
    
    Returns:
        List of matching GFF entries
    """
    with get_starbase_session() as session:
        query = (
            session.query(Gff, Ships, ShipAccessions)
            .join(Ships, Gff.ship_id == Ships.id)
            .outerjoin(ShipAccessions, ShipAccessions.ship_id == Ships.id)
        )
        if ship_ids:
            query = query.filter(Gff.ship_id.in_(ship_ids))
        query = query.filter(Gff.attributes.like(f'%Target_ID={target_id}%'))
        query = query.order_by(Ships.id, Gff.start)
        
        results = []
        for gff, ship, ship_acc in query.all():
            results.append({
                'id': gff.id,
                'start': gff.start,
                'end': gff.end,
                'strand': gff.strand,
                'type': gff.type,
                'attributes': gff.attributes,
                'ship_id': gff.ship_id,
                'ship_accession': ship_acc.ship_accession_display if ship_acc else '',
            })
        
        return results


def get_genes_in_region(ship_id, start, end):
    """
    Get all genes in a specific genomic region.
    
    Args:
        ship_id: Ship ID
        start: Start position
        end: End position
    
    Returns:
        List of GFF entries in the region
    """
    with get_starbase_session() as session:
        query = (
            session.query(Gff)
            .filter(Gff.ship_id == ship_id)
            .filter(
                ((Gff.start >= start) & (Gff.start <= end)) |
                ((Gff.end >= start) & (Gff.end <= end)) |
                ((Gff.start <= start) & (Gff.end >= end))
            )
            .order_by(Gff.start)
        )
        
        results = []
        for gff in query.all():
            results.append({
                'id': gff.id,
                'start': gff.start,
                'end': gff.end,
                'strand': gff.strand,
                'type': gff.type,
                'source': gff.source,
                'attributes': gff.attributes,
            })
        
        return results
