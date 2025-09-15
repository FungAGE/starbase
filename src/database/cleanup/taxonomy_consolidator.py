#!/usr/bin/env python3
"""
Taxonomy Consolidator

This module provides a comprehensive, consolidated approach to taxonomy table
management, including data cleaning, ID resolution, hierarchy backfill, and
consistency validation.
"""

import re
import logging
from typing import Dict, List, Optional, Set, Tuple, Any
from sqlalchemy import text, func
from sqlalchemy.orm import sessionmaker

from ..models.schema import Taxonomy, Genome, JoinedShips
from ..sql_manager import StarbaseSession

logger = logging.getLogger(__name__)


class TaxonomyConsolidator:
    """
    Comprehensive taxonomy table consolidation and management.
    
    This class provides a unified approach to:
    1. Data cleaning and integrity checks
    2. Taxonomy ID resolution from multiple sources
    3. Taxonomy hierarchy backfill
    4. Consistency validation
    """
    
    def __init__(self, session=None):
        """
        Initialize the TaxonomyConsolidator.
        
        Args:
            session: Database session (optional, creates new if not provided)
        """
        self.session = session or StarbaseSession()
        self._owns_session = session is None
        
        # Taxonomy hierarchy fields in order from highest to lowest
        self.hierarchy_fields = [
            'superkingdom', 'clade', 'kingdom', 'subkingdom', 'phylum', 
            'subphylum', 'class_', 'subclass', 'order', 'suborder', 
            'family', 'genus', 'species', 'section', 'species_group', 
            'subgenus', 'strain'
        ]
        
        # Fields that should be preserved (not backfilled)
        self.preserved_fields = {
            'name', 'taxID', 'genus', 'species', 'section', 
            'species_group', 'subgenus', 'strain'
        }
        
        # Fields that can be backfilled
        self.backfillable_fields = [
            'superkingdom', 'clade', 'kingdom', 'subkingdom', 'phylum',
            'subphylum', 'class_', 'subclass', 'order', 'suborder', 'family'
        ]
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._owns_session:
            self.session.close()
    
    def run_full_consolidation(self, 
                              clean_data: bool = True,
                              resolve_ids: bool = True, 
                              backfill_hierarchy: bool = True,
                              validate_consistency: bool = True,
                              dry_run: bool = True,
                              names_dmp_path: str = None,
                              ncbi_email: str = None,
                              ncbi_api_key: str = None,
                              include_synonyms: bool = False) -> Dict:
        """
        Run complete taxonomy consolidation pipeline.
        Args:
            clean_data: Whether to clean taxonomy data                                                                                                                                                                                                                                                                                  
            resolve_ids: Whether to resolve missing taxonomy IDs
            backfill_hierarchy: Whether to backfill missing hierarchy fields
            validate_consistency: Whether to validate consistency
            dry_run: If True, only analyze and report what would be done
            
        Returns:                            
            Dict: Comprehensive report of all operations performed
        """
        logger.info("Starting comprehensive taxonomy consolidation...")
        
        full_report = {
            'cleaning': {},
            'id_resolution': {},
            'hierarchy_backfill': {},
            'consistency_validation': {},
            'summary': {}
        }
        
        try:
            # Phase 1: Data cleaning
            if clean_data:
                logger.info("Phase 1: Cleaning taxonomy data...")
                full_report['cleaning'] = self.clean_taxonomy_data(dry_run=dry_run)
            
            # Phase 2: Taxonomy ID resolution
            if resolve_ids:
                logger.info("Phase 2: Resolving taxonomy IDs...")
                full_report['id_resolution'] = self.resolve_taxonomy_ids(
                    dry_run=dry_run,
                    names_dmp_path=names_dmp_path,
                    ncbi_email=ncbi_email,
                    ncbi_api_key=ncbi_api_key,
                    include_synonyms=include_synonyms
                )
            
            # Phase 3: Hierarchy backfill
            if backfill_hierarchy:
                logger.info("Phase 3: Backfilling taxonomy hierarchy...")
                full_report['hierarchy_backfill'] = self.backfill_taxonomy_hierarchy(dry_run=dry_run)
            
            # Phase 4: Consistency validation
            if validate_consistency:
                logger.info("Phase 4: Validating taxonomy consistency...")
                full_report['consistency_validation'] = self.validate_taxonomy_consistency()
            
            # Generate summary
            full_report['summary'] = self._generate_consolidation_summary(full_report)
            
            if not dry_run:
                self.session.commit()
                logger.info("Taxonomy consolidation completed successfully")
            else:
                logger.info("Taxonomy consolidation analysis completed (dry run)")
            
            return full_report
            
        except Exception as e:
            logger.error(f"Error during taxonomy consolidation: {str(e)}")
            self.session.rollback()
            raise
    
    def clean_taxonomy_data(self, dry_run: bool = True) -> Dict:
        """
        Phase 1: Clean taxonomy data and fix integrity issues.
        
        Tasks:
        1. Remove trailing/leading whitespace, normalize line endings
        2. Fix duplicated information in species/strain fields
        3. Clean up name field duplications
        4. Validate taxonomic hierarchy consistency
        5. Flag data quality issues
        
        Args:
            dry_run: If True, only analyze and report what would be done
            
        Returns:
            Dict: Report of cleaning operations performed
        """
        logger.info("Cleaning taxonomy data...")
        
        report = {
            'whitespace_fixes': [],
            'species_cleanup': [],
            'strain_cleanup': [],
            'name_cleanup': [],
            'hierarchy_issues': [],
            'summary': {}
        }
        
        try:
            # Get all taxonomy entries
            taxonomies = self.session.query(Taxonomy).all()
            logger.info(f"Processing {len(taxonomies)} taxonomy entries for cleaning")
            
            whitespace_fixes = 0
            species_cleanup = 0
            strain_cleanup = 0
            name_cleanup = 0
            hierarchy_issues = 0
            
            for tax in taxonomies:
                changes_made = False
                
                # 1. Convert empty strings to NULL and fix whitespace issues
                for field in self.hierarchy_fields + ['name', 'taxID']:
                    old_value = getattr(tax, field)
                    if old_value == "":
                        # Convert empty string to None (NULL)
                        if not dry_run:
                            setattr(tax, field, None)
                        report['whitespace_fixes'].append({
                            'taxonomy_id': tax.id,
                            'field': field,
                            'old_value': '""',
                            'new_value': 'NULL'
                        })
                        changes_made = True
                        whitespace_fixes += 1
                    elif old_value:
                        # Normalize whitespace and line endings
                        new_value = self._normalize_whitespace(old_value)
                        if new_value != old_value:
                            if not dry_run:
                                setattr(tax, field, new_value)
                            report['whitespace_fixes'].append({
                                'taxonomy_id': tax.id,
                                'field': field,
                                'old_value': old_value,
                                'new_value': new_value
                            })
                            changes_made = True
                            whitespace_fixes += 1
                
                # 2. Clean species field (remove genus if duplicated)
                if tax.species and tax.genus:
                    original_species = tax.species
                    cleaned_species = self._clean_species_field(tax.species, tax.genus)
                    if cleaned_species != original_species:
                        if not dry_run:
                            tax.species = cleaned_species
                        report['species_cleanup'].append({
                            'taxonomy_id': tax.id,
                            'old_species': original_species,
                            'new_species': cleaned_species,
                            'genus': tax.genus
                        })
                        changes_made = True
                        species_cleanup += 1
                
                # 3. Clean strain field (remove genus/species matches)
                if tax.strain and (tax.genus or tax.species):
                    original_strain = tax.strain
                    cleaned_strain = self._clean_strain_field(tax.strain, tax.genus, tax.species)
                    if cleaned_strain != original_strain:
                        if not dry_run:
                            tax.strain = cleaned_strain
                        report['strain_cleanup'].append({
                            'taxonomy_id': tax.id,
                            'old_strain': original_strain,
                            'new_strain': cleaned_strain,
                            'genus': tax.genus,
                            'species': tax.species
                        })
                        changes_made = True
                        strain_cleanup += 1
                
                # 4. Set name field to concatenation of genus + species + strain
                expected_name = self._build_expected_name(tax.genus, tax.species, tax.strain)
                if expected_name and tax.name != expected_name:
                    original_name = tax.name
                    if not dry_run:
                        tax.name = expected_name
                    report['name_cleanup'].append({
                        'taxonomy_id': tax.id,
                        'old_name': original_name,
                        'new_name': expected_name,
                        'genus': tax.genus,
                        'species': tax.species,
                        'strain': tax.strain
                    })
                    changes_made = True
                    name_cleanup += 1
                
                # 5. Check for hierarchy consistency issues
                hierarchy_issue = self._check_hierarchy_consistency(tax)
                if hierarchy_issue:
                    report['hierarchy_issues'].append({
                        'taxonomy_id': tax.id,
                        'name': tax.name,
                        'issue': hierarchy_issue
                    })
                    hierarchy_issues += 1
                
                if changes_made and not dry_run:
                    self.session.add(tax)
            
            report['summary'] = {
                'total_processed': len(taxonomies),
                'whitespace_fixes': whitespace_fixes,
                'species_cleanup': species_cleanup,
                'strain_cleanup': strain_cleanup,
                'name_cleanup': name_cleanup,
                'hierarchy_issues': hierarchy_issues
            }
            
            logger.info(f"Cleaning completed: {whitespace_fixes} whitespace fixes, "
                       f"{species_cleanup} species cleanups, {strain_cleanup} strain cleanups, "
                       f"{name_cleanup} name cleanups, {hierarchy_issues} hierarchy issues")
            
            # Commit changes if not dry run
            if not dry_run and (whitespace_fixes > 0 or species_cleanup > 0 or 
                               strain_cleanup > 0 or name_cleanup > 0):
                self.session.commit()
                logger.info("Changes committed to database")
            
            return report
            
        except Exception as e:
            logger.error(f"Error during taxonomy data cleaning: {str(e)}")
            raise
    
    def resolve_taxonomy_ids(self, dry_run: bool = True, names_dmp_path: str = None, 
                           ncbi_email: str = None, ncbi_api_key: str = None, 
                           include_synonyms: bool = False) -> Dict:
        """
        Phase 2: Resolve missing taxonomy IDs from multiple sources.
        
        Priority order:
        1. Existing taxIDs from other tables (joined_ships, etc.)
        2. NCBI taxdump names.dmp lookup
        3. NCBI API lookup (with rate limiting)
        4. Flag unresolved entries for manual review
        
        Args:
            dry_run: If True, only analyze and report what would be done
            
        Returns:
            Dict: Report of ID resolution operations
        """
        logger.info("Resolving taxonomy IDs...")
        
        # Initialize taxdump and NCBI credentials
        self._taxdump_names = None
        self._ncbi_email = ncbi_email
        self._ncbi_api_key = ncbi_api_key
        
        # Load taxdump data if path provided
        if names_dmp_path:
            try:
                from .utils.database_cleanup import _parse_taxdump_names
                self._taxdump_names = _parse_taxdump_names(names_dmp_path, include_synonyms)
                logger.info(f"Loaded {len(self._taxdump_names)} names from taxdump")
            except Exception as e:
                logger.warning(f"Failed to load taxdump: {str(e)}")
                self._taxdump_names = None
        
        report = {
            'resolved_from_existing': [],
            'resolved_from_taxdump': [],
            'resolved_from_ncbi': [],
            'unresolved': [],
            'summary': {}
        }
        
        try:
            # Find taxonomies without taxID
            missing_taxid = self.session.query(Taxonomy).filter(
                (Taxonomy.taxID.is_(None)) | (Taxonomy.taxID == '')
            ).all()
            
            logger.info(f"Found {len(missing_taxid)} taxonomies without taxID")
            
            resolved_count = 0
            
            for tax in missing_taxid:
                taxid = None
                source = None
                
                # 1. Try to find existing taxID from other tables
                taxid = self._find_existing_taxid(tax)
                if taxid:
                    source = 'existing'
                    report['resolved_from_existing'].append({
                        'taxonomy_id': tax.id,
                        'name': tax.name,
                        'taxID': taxid,
                        'source': source
                    })
                
                # 2. Try taxdump lookup (if available)
                if not taxid:
                    taxid = self._lookup_taxid_from_taxdump(tax)
                    if taxid:
                        source = 'taxdump'
                        report['resolved_from_taxdump'].append({
                            'taxonomy_id': tax.id,
                            'name': tax.name,
                            'taxID': taxid,
                            'source': source
                        })
                
                # 3. Try NCBI API lookup (if available)
                if not taxid:
                    taxid = self._lookup_taxid_from_ncbi(tax)
                    if taxid:
                        source = 'ncbi'
                        report['resolved_from_ncbi'].append({
                            'taxonomy_id': tax.id,
                            'name': tax.name,
                            'taxID': taxid,
                            'source': source
                        })
                
                # Update taxonomy with resolved taxID
                if taxid:
                    if not dry_run:
                        tax.taxID = str(taxid)
                        self.session.add(tax)
                    resolved_count += 1
                else:
                    report['unresolved'].append({
                        'taxonomy_id': tax.id,
                        'name': tax.name,
                        'genus': tax.genus,
                        'species': tax.species,
                        'reason': 'No taxID found in any source'
                    })
            
            report['summary'] = {
                'total_missing_taxid': len(missing_taxid),
                'resolved': resolved_count,
                'unresolved': len(report['unresolved']),
                'resolved_from_existing': len(report['resolved_from_existing']),
                'resolved_from_taxdump': len(report['resolved_from_taxdump']),
                'resolved_from_ncbi': len(report['resolved_from_ncbi'])
            }
            
            logger.info(f"ID resolution completed: {resolved_count} resolved, "
                       f"{len(report['unresolved'])} unresolved")
            
            return report
            
        except Exception as e:
            logger.error(f"Error during taxonomy ID resolution: {str(e)}")
            raise
    
    def backfill_taxonomy_hierarchy(self, dry_run: bool = True) -> Dict:
        """
        Phase 3: Backfill missing taxonomy hierarchy fields.
        
        Args:
            dry_run: If True, only analyze and report what would be done
            
        Returns:
            Dict: Report of hierarchy backfill operations
        """
        logger.info("Backfilling taxonomy hierarchy...")
        
        report = {
            'fields_backfilled': [],
            'consistency_issues': [],
            'summary': {}
        }
        
        try:
            # Get taxonomies with missing hierarchy fields
            taxonomies = self.session.query(Taxonomy).all()
            
            backfilled_count = 0
            consistency_issues = 0
            
            for tax in taxonomies:
                changes_made = False
                
                # Try to backfill missing fields from existing data
                for field in self.backfillable_fields:
                    current_value = getattr(tax, field)
                    if not current_value or current_value.strip() == '':
                        # Try to infer from other entries with same family/genus
                        inferred_value = self._infer_hierarchy_field(tax, field)
                        if inferred_value:
                            if not dry_run:
                                setattr(tax, field, inferred_value)
                            report['fields_backfilled'].append({
                                'taxonomy_id': tax.id,
                                'field': field,
                                'inferred_value': inferred_value,
                                'source': 'family_consistency'
                            })
                            changes_made = True
                            backfilled_count += 1
                
                # Check for consistency issues
                consistency_issue = self._check_field_consistency(tax)
                if consistency_issue:
                    report['consistency_issues'].append({
                        'taxonomy_id': tax.id,
                        'name': tax.name,
                        'issue': consistency_issue
                    })
                    consistency_issues += 1
                
                if changes_made and not dry_run:
                    self.session.add(tax)
            
            # Commit changes if not in dry-run mode
            if not dry_run and backfilled_count > 0:
                self.session.commit()
                logger.info(f"Committed {backfilled_count} backfilled fields to database")
            
            report['summary'] = {
                'total_processed': len(taxonomies),
                'fields_backfilled': backfilled_count,
                'consistency_issues': consistency_issues
            }
            
            logger.info(f"Hierarchy backfill completed: {backfilled_count} fields backfilled, "
                       f"{consistency_issues} consistency issues found")
            
            return report
            
        except Exception as e:
            logger.error(f"Error during taxonomy hierarchy backfill: {str(e)}")
            raise
    
    def validate_taxonomy_consistency(self) -> Dict:
        """
        Phase 4: Comprehensive taxonomy consistency validation.
        
        Returns:
            Dict: Report of consistency validation results
        """
        logger.info("Validating taxonomy consistency...")
        
        report = {
            'family_inconsistencies': [],
            'hierarchy_violations': [],
            'orphaned_entries': [],
            'duplicate_entries': [],
            'data_quality_issues': [],
            'summary': {}
        }
        
        try:
            # 1. Check family consistency
            family_issues = self._check_family_consistency()
            report['family_inconsistencies'] = family_issues
            
            # 2. Check hierarchy violations
            hierarchy_violations = self._check_hierarchy_violations()
            report['hierarchy_violations'] = hierarchy_violations
            
            # 3. Check for orphaned entries
            orphaned = self._check_orphaned_entries()
            report['orphaned_entries'] = orphaned
            
            # 4. Check for duplicates
            duplicates = self._check_duplicate_entries()
            report['duplicate_entries'] = duplicates
            
            # 5. Data quality assessment
            quality_issues = self._assess_data_quality()
            report['data_quality_issues'] = quality_issues
            
            report['summary'] = {
                'family_inconsistencies': len(family_issues),
                'hierarchy_violations': len(hierarchy_violations),
                'orphaned_entries': len(orphaned),
                'duplicate_entries': len(duplicates),
                'data_quality_issues': len(quality_issues)
            }
            
            logger.info(f"Consistency validation completed: "
                       f"{len(family_issues)} family inconsistencies, "
                       f"{len(hierarchy_violations)} hierarchy violations, "
                       f"{len(orphaned)} orphaned entries, "
                       f"{len(duplicates)} duplicate entries, "
                       f"{len(quality_issues)} data quality issues")
            
            return report
            
        except Exception as e:
            logger.error(f"Error during taxonomy consistency validation: {str(e)}")
            raise
    
    # Helper methods
    
    def _normalize_whitespace(self, text: str) -> str:
        """Normalize whitespace and line endings in text."""
        if not text:
            return text
        
        # Replace various line ending types with single space
        text = re.sub(r'[\r\n]+', ' ', text)
        # Collapse multiple spaces
        text = re.sub(r'\s+', ' ', text)
        # Strip leading/trailing whitespace
        return text.strip()
    
    def _clean_species_field(self, species: str, genus: str) -> str:
        """Remove genus from species field if duplicated."""
        if not species or not genus:
            return species
        
        # Normalize both strings for comparison
        species_lower = species.lower().strip()
        genus_lower = genus.lower().strip()
        
        # Case 1: Species starts with genus (e.g., "Alternaria alstroemeriae" -> "alstroemeriae")
        if species_lower.startswith(genus_lower):
            # Remove genus and any following space
            cleaned = species[len(genus):].strip()
            return cleaned if cleaned else species
        
        # Case 2: Species contains genus as a word (e.g., "Alternaria alstroemeriae" -> "alstroemeriae")
        # This handles cases where there might be extra spaces or formatting
        genus_pattern = re.compile(r'\b' + re.escape(genus_lower) + r'\b', re.IGNORECASE)
        if genus_pattern.search(species):
            # Remove the genus word and clean up extra spaces
            cleaned = genus_pattern.sub('', species).strip()
            # Clean up any double spaces
            cleaned = re.sub(r'\s+', ' ', cleaned).strip()
            return cleaned if cleaned else species
        
        return species
    
    def _clean_strain_field(self, strain: str, genus: str, species: str) -> str:
        """Remove genus/species matches from strain field."""
        if not strain:
            return strain
        
        cleaned = strain
        
        # Remove full genus matches
        if genus:
            # Remove exact genus matches (case insensitive)
            pattern = re.compile(r'\b' + re.escape(genus) + r'\b', re.IGNORECASE)
            cleaned = pattern.sub('', cleaned).strip()
        
        # Remove full species matches
        if species:
            pattern = re.compile(r'\b' + re.escape(species) + r'\b', re.IGNORECASE)
            cleaned = pattern.sub('', cleaned).strip()
        
        # Clean up extra spaces
        cleaned = re.sub(r'\s+', ' ', cleaned).strip()
        
        return cleaned if cleaned else strain
    
    def _clean_name_field(self, name: str, genus: str, species: str) -> str:
        """Remove duplicated genus+species from name field."""
        if not name or not genus or not species:
            return name
        
        # Create expected genus+species combination
        expected_combination = f"{genus} {species}"
        
        # Check for duplicated genus+species in name
        if expected_combination.lower() in name.lower():
            # Find and remove the first occurrence
            pattern = re.compile(re.escape(expected_combination), re.IGNORECASE)
            cleaned = pattern.sub('', name, count=1).strip()
            # Clean up extra spaces
            cleaned = re.sub(r'\s+', ' ', cleaned).strip()
            return cleaned
        
        return name
    
    def _build_expected_name(self, genus: str, species: str, strain: str) -> str:
        """Build expected name field from genus + species + strain."""
        parts = []
        
        if genus and genus.strip():
            parts.append(genus.strip())
        
        if species and species.strip():
            parts.append(species.strip())
        
        if strain and strain.strip():
            parts.append(strain.strip())
        
        if parts:
            return ' '.join(parts)
        
        return None
    
    def _check_hierarchy_consistency(self, tax: Taxonomy) -> Optional[str]:
        """Check for hierarchy consistency issues in a taxonomy entry."""
        issues = []
        
        # Check if genus is present but family is missing
        if tax.genus and not tax.family:
            issues.append("Genus present but family missing")
        
        # Check if species is present but genus is missing
        if tax.species and not tax.genus:
            issues.append("Species present but genus missing")
        
        # Check if family is present but order is missing
        if tax.family and not tax.order:
            issues.append("Family present but order missing")
        
        return "; ".join(issues) if issues else None
    
    def _find_existing_taxid(self, tax: Taxonomy) -> Optional[str]:
        """Find existing taxID from other tables."""
        # Look in joined_ships table
        joined_entry = self.session.query(JoinedShips).filter(
            JoinedShips.tax_id == tax.id
        ).first()
        
        if joined_entry:
            # Try to find a related taxonomy with taxID
            related_tax = self.session.query(Taxonomy).filter(
                (Taxonomy.genus == tax.genus) & 
                (Taxonomy.species == tax.species) &
                (Taxonomy.taxID.isnot(None)) &
                (Taxonomy.taxID != '')
            ).first()
            
            if related_tax:
                return related_tax.taxID
        
        return None
    
    def _lookup_taxid_from_taxdump(self, tax: Taxonomy) -> Optional[str]:
        """Lookup taxID from taxdump using existing functionality."""
        if not hasattr(self, '_taxdump_names') or self._taxdump_names is None:
            return None
        
        # Try to match by name, genus+species, or strain
        search_terms = []
        
        if tax.name and tax.name.strip():
            search_terms.append(tax.name.strip())
        
        if tax.genus and tax.species:
            search_terms.append(f"{tax.genus} {tax.species}")
        
        if tax.strain and tax.strain.strip():
            search_terms.append(tax.strain.strip())
        
        for term in search_terms:
            # Try exact match
            if term in self._taxdump_names:
                return self._taxdump_names[term]
            
            # Try case-insensitive match
            term_lower = term.lower()
            for name, taxid in self._taxdump_names.items():
                if name.lower() == term_lower:
                    return taxid
        
        return None
    
    def _lookup_taxid_from_ncbi(self, tax: Taxonomy) -> Optional[str]:
        """Lookup taxID from NCBI API using existing functionality."""
        if not hasattr(self, '_ncbi_email') or not self._ncbi_email:
            return None
        
        # Import the existing NCBI function
        from .utils.database_cleanup import _ncbi_taxonomy_search
        
        # Try to match by name, genus+species, or strain
        search_terms = []
        
        if tax.name and tax.name.strip():
            search_terms.append(tax.name.strip())
        
        if tax.genus and tax.species:
            search_terms.append(f"{tax.genus} {tax.species}")
        
        for term in search_terms:
            try:
                result = _ncbi_taxonomy_search(
                    scientific_name=term,
                    email=self._ncbi_email,
                    api_key=getattr(self, '_ncbi_api_key', None)
                )
                if result and 'taxid' in result:
                    return result['taxid']
            except Exception as e:
                logger.warning(f"NCBI lookup failed for '{term}': {str(e)}")
                continue
        
        return None
    
    def _infer_hierarchy_field(self, tax: Taxonomy, field: str) -> Optional[str]:
        """Infer missing hierarchy field from other entries with same family/genus."""
        if not tax.family and not tax.genus:
            return None
        
        # Strategy 1: Find entries with same family (most specific)
        if tax.family:
            matching_entries = self.session.query(Taxonomy).filter(
                Taxonomy.id != tax.id,
                Taxonomy.family == tax.family,
                getattr(Taxonomy, field).isnot(None),
                getattr(Taxonomy, field) != ''
            ).all()
            
            if matching_entries:
                # Return the most common value
                values = [getattr(entry, field) for entry in matching_entries]
                return max(set(values), key=values.count)
        
        # Strategy 2: Find entries with same genus (less specific but still useful)
        if tax.genus:
            matching_entries = self.session.query(Taxonomy).filter(
                Taxonomy.id != tax.id,
                Taxonomy.genus == tax.genus,
                getattr(Taxonomy, field).isnot(None),
                getattr(Taxonomy, field) != ''
            ).all()
            
            if matching_entries:
                # Return the most common value
                values = [getattr(entry, field) for entry in matching_entries]
                return max(set(values), key=values.count)
        
        return None
    
    def _check_field_consistency(self, tax: Taxonomy) -> Optional[str]:
        """Check for field consistency issues."""
        issues = []
        
        # Check if name matches genus+species
        if tax.name and tax.genus and tax.species:
            expected_name = f"{tax.genus} {tax.species}"
            if tax.name.lower() != expected_name.lower():
                issues.append(f"Name '{tax.name}' doesn't match genus+species '{expected_name}'")
        
        return "; ".join(issues) if issues else None
    
    def _check_family_consistency(self) -> List[Dict]:
        """Check for family consistency issues."""
        issues = []
        
        # Group by family and check for inconsistencies
        families = self.session.query(Taxonomy.family).filter(
            Taxonomy.family.isnot(None),
            Taxonomy.family != ''
        ).distinct().all()
        
        for (family,) in families:
            family_entries = self.session.query(Taxonomy).filter(
                Taxonomy.family == family
            ).all()
            
            # Check for inconsistent upstream taxonomy
            for field in ['superkingdom', 'kingdom', 'phylum', 'class_', 'order']:
                values = [getattr(entry, field) for entry in family_entries 
                         if getattr(entry, field) and getattr(entry, field).strip()]
                
                if len(set(values)) > 1:  # Multiple different values
                    issues.append({
                        'family': family,
                        'field': field,
                        'values': list(set(values)),
                        'count': len(family_entries)
                    })
        
        return issues
    
    def _check_hierarchy_violations(self) -> List[Dict]:
        """Check for taxonomic hierarchy violations."""
        violations = []
        
        # This would implement more sophisticated hierarchy validation
        # For now, return empty list
        return violations
    
    def _check_orphaned_entries(self) -> List[Dict]:
        """Check for orphaned taxonomy entries."""
        orphaned = []
        
        # Find taxonomies with no associated genomes or joined_ships
        orphaned_entries = self.session.query(Taxonomy).outerjoin(
            Genome, Taxonomy.id == Genome.taxonomy_id
        ).outerjoin(
            JoinedShips, Taxonomy.id == JoinedShips.tax_id
        ).filter(
            Genome.id.is_(None),
            JoinedShips.id.is_(None)
        ).all()
        
        for entry in orphaned_entries:
            orphaned.append({
                'taxonomy_id': entry.id,
                'name': entry.name,
                'taxID': entry.taxID
            })
        
        return orphaned
    
    def _check_duplicate_entries(self) -> List[Dict]:
        """Check for duplicate taxonomy entries based on name, genus, species, and strain."""
        duplicates = []
        
        # Find entries with same combination of name, genus, species, and strain
        # This is more reliable than taxID since taxIDs can be upstream NCBI IDs
        duplicate_combinations = self.session.query(
            Taxonomy.name, Taxonomy.genus, Taxonomy.species, Taxonomy.strain
        ).filter(
            # Only check entries that have at least name or genus+species
            (Taxonomy.name.isnot(None) & (Taxonomy.name != '')) |
            ((Taxonomy.genus.isnot(None) & (Taxonomy.genus != '')) & 
             (Taxonomy.species.isnot(None) & (Taxonomy.species != '')))
        ).group_by(
            Taxonomy.name, Taxonomy.genus, Taxonomy.species, Taxonomy.strain
        ).having(
            func.count(Taxonomy.id) > 1
        ).all()
        
        for (name, genus, species, strain) in duplicate_combinations:
            entries = self.session.query(Taxonomy).filter(
                Taxonomy.name == name,
                Taxonomy.genus == genus,
                Taxonomy.species == species,
                Taxonomy.strain == strain
            ).all()
            
            duplicates.append({
                'name': name,
                'genus': genus,
                'species': species,
                'strain': strain,
                'count': len(entries),
                'entries': [{'id': e.id, 'name': e.name, 'taxID': e.taxID} for e in entries]
            })
        
        return duplicates
    
    def consolidate_duplicate_entries(self, dry_run: bool = True) -> Dict:
        """
        Consolidate duplicate taxonomy entries and update foreign keys.
        
        Args:
            dry_run: If True, only report what would be done without making changes
            
        Returns:
            Dict: Report of consolidation actions
        """
        logger.info("Consolidating duplicate taxonomy entries...")
        
        report = {
            'duplicates_found': [],
            'consolidations_performed': [],
            'foreign_key_updates': [],
            'summary': {}
        }
        
        try:
            # Find all duplicates
            duplicates = self._check_duplicate_entries()
            report['duplicates_found'] = duplicates
            
            consolidations = 0
            fk_updates = 0
            
            for duplicate_group in duplicates:
                entries = self.session.query(Taxonomy).filter(
                    Taxonomy.name == duplicate_group['name'],
                    Taxonomy.genus == duplicate_group['genus'],
                    Taxonomy.species == duplicate_group['species'],
                    Taxonomy.strain == duplicate_group['strain']
                ).all()
                
                if len(entries) <= 1:
                    continue
                
                # Sort entries to determine which one to keep
                # Priority: 1) Has taxID, 2) Has more complete data, 3) Lower ID (older)
                def entry_priority(entry):
                    score = 0
                    if entry.taxID and entry.taxID.strip():
                        score += 1000
                    # Count non-empty fields
                    non_empty_fields = sum(1 for field in [
                        entry.name, entry.genus, entry.species, entry.strain,
                        entry.family, entry.order, entry.phylum, entry.kingdom
                    ] if field and field.strip())
                    score += non_empty_fields * 10
                    # Lower ID is better (older entry)
                    score += (10000 - entry.id)
                    return score
                
                entries.sort(key=entry_priority, reverse=True)
                keep_entry = entries[0]
                remove_entries = entries[1:]
                
                consolidation_info = {
                    'name': duplicate_group['name'],
                    'genus': duplicate_group['genus'],
                    'species': duplicate_group['species'],
                    'strain': duplicate_group['strain'],
                    'keep_entry': {'id': keep_entry.id, 'name': keep_entry.name, 'taxID': keep_entry.taxID},
                    'remove_entries': [{'id': e.id, 'name': e.name, 'taxID': e.taxID} for e in remove_entries]
                }
                
                if not dry_run:
                    # Update foreign keys in JoinedShips
                    for remove_entry in remove_entries:
                        joined_ships_updates = self.session.query(JoinedShips).filter(
                            JoinedShips.tax_id == remove_entry.id
                        ).update({JoinedShips.tax_id: keep_entry.id})
                        
                        if joined_ships_updates > 0:
                            report['foreign_key_updates'].append({
                                'table': 'joined_ships',
                                'old_taxonomy_id': remove_entry.id,
                                'new_taxonomy_id': keep_entry.id,
                                'records_updated': joined_ships_updates
                            })
                            fk_updates += joined_ships_updates
                    
                    # Update foreign keys in Genomes
                    for remove_entry in remove_entries:
                        genome_updates = self.session.query(Genome).filter(
                            Genome.taxonomy_id == remove_entry.id
                        ).update({Genome.taxonomy_id: keep_entry.id})
                        
                        if genome_updates > 0:
                            report['foreign_key_updates'].append({
                                'table': 'genomes',
                                'old_taxonomy_id': remove_entry.id,
                                'new_taxonomy_id': keep_entry.id,
                                'records_updated': genome_updates
                            })
                            fk_updates += genome_updates
                    
                    # Remove duplicate entries
                    for remove_entry in remove_entries:
                        self.session.delete(remove_entry)
                    
                    consolidations += len(remove_entries)
                
                report['consolidations_performed'].append(consolidation_info)
            
            report['summary'] = {
                'duplicates_found': len(duplicates),
                'consolidations_performed': consolidations,
                'foreign_key_updates': fk_updates
            }
            
            logger.info(f"Consolidation completed: {len(duplicates)} duplicate groups found, "
                       f"{consolidations} entries consolidated, {fk_updates} foreign keys updated")
            
            # Commit changes if not dry run
            if not dry_run and (consolidations > 0 or fk_updates > 0):
                self.session.commit()
                logger.info("Duplicate consolidation changes committed to database")
            
            return report
            
        except Exception as e:
            logger.error(f"Error during duplicate consolidation: {str(e)}")
            raise
    
    def _assess_data_quality(self) -> List[Dict]:
        """Assess overall data quality."""
        issues = []
        
        # Check for entries with missing essential fields
        missing_essential = self.session.query(Taxonomy).filter(
            (Taxonomy.name.is_(None)) | (Taxonomy.name == '')
        ).count()
        
        if missing_essential > 0:
            issues.append({
                'type': 'missing_essential_fields',
                'count': missing_essential,
                'description': 'Entries missing name field'
            })
        
        return issues
    
    def _generate_consolidation_summary(self, full_report: Dict) -> Dict:
        """Generate summary of consolidation operations."""
        summary = {
            'total_operations': 0,
            'successful_operations': 0,
            'issues_found': 0,
            'recommendations': []
        }
        
        # Count operations and issues
        for phase, report in full_report.items():
            if isinstance(report, dict) and 'summary' in report:
                phase_summary = report['summary']
                if isinstance(phase_summary, dict):
                    summary['total_operations'] += phase_summary.get('total_processed', 0)
                    summary['successful_operations'] += sum(
                        v for k, v in phase_summary.items() 
                        if k.endswith('_fixes') or k.endswith('_cleanup') or k == 'resolved'
                    )
                    summary['issues_found'] += sum(
                        v for k, v in phase_summary.items() 
                        if k.endswith('_issues') or k == 'unresolved'
                    )
        
        # Generate recommendations
        if summary['issues_found'] > 0:
            summary['recommendations'].append(
                f"Found {summary['issues_found']} issues that may require manual review"
            )
        
        if summary['successful_operations'] > 0:
            summary['recommendations'].append(
                f"Successfully processed {summary['successful_operations']} operations"
            )
        
        return summary
