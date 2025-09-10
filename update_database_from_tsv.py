#!/usr/bin/env python3
"""
Script to update SQLite database from edited TSV file.
This script safely updates the database based on changes in the TSV export.
"""

import csv
import sys
import os
from pathlib import Path
from typing import Dict, List, Optional, Any

# Add src to path so we can import database modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from database.sql_engine import get_starbase_session
from database.models.schema import (
    Accessions, JoinedShips, Taxonomy, FamilyNames, 
    NavisNames, HaplotypeNames, Genome, Ships, Base
)
from sqlalchemy import and_
from config.logging import get_logger

logger = get_logger(__name__)


class DatabaseUpdater:
    """Handle database updates from TSV data"""
    
    def __init__(self, tsv_file_path: str):
        self.tsv_file_path = tsv_file_path
        self.stats = {
            'processed': 0,
            'updated': 0,
            'errors': 0,
            'skipped': 0
        }
    
    def read_tsv_data(self) -> List[Dict[str, Any]]:
        """Read TSV file and return list of dictionaries"""
        data = []
        try:
            with open(self.tsv_file_path, 'r', encoding='utf-8') as file:
                reader = csv.DictReader(file, delimiter='\t')
                for row in reader:
                    # Clean up empty strings and convert to None where appropriate
                    cleaned_row = {}
                    for key, value in row.items():
                        if value == '' or value == '.':
                            cleaned_row[key] = None
                        else:
                            cleaned_row[key] = value.strip() if isinstance(value, str) else value
                    data.append(cleaned_row)
            logger.info(f"Successfully read {len(data)} rows from TSV file")
            return data
        except Exception as e:
            logger.error(f"Error reading TSV file: {e}")
            raise
    
    def update_accession(self, session, row: Dict[str, Any]) -> Optional[Accessions]:
        """Update or create accession record"""
        accession_tag = row.get('accession_tag')
        if not accession_tag:
            return None
            
        accession = session.query(Accessions).filter(
            Accessions.accession_tag == accession_tag
        ).first()
        
        if accession:
            # Update existing accession if needed
            # Most accession data seems to be stable, but we can add updates here if needed
            pass
        else:
            # Create new accession with default version_tag = 1
            accession = Accessions(
                accession_tag=accession_tag,
                ship_name=row.get('starshipID'),
                version_tag=1  # Set default version_tag to satisfy NOT NULL constraint
            )
            session.add(accession)
            session.flush()  # Get the ID
            logger.info(f"Created new accession: {accession_tag}")
        
        return accession
    
    def update_taxonomy(self, session, row: Dict[str, Any]) -> Optional[Taxonomy]:
        """Update or create taxonomy record"""
        taxid = row.get('taxID')
        if not taxid:
            return None
            
        try:
            taxid = int(taxid)
        except (ValueError, TypeError):
            logger.warning(f"Invalid taxID: {taxid}")
            return None
        
        taxonomy = session.query(Taxonomy).filter(
            Taxonomy.taxID == str(taxid)
        ).first()
        
        if taxonomy:
            # Update existing taxonomy
            if row.get('order'):
                taxonomy.order = row['order']
            if row.get('family'):
                taxonomy.family = row['family']
            if row.get('genus'):
                taxonomy.genus = row['genus']
            if row.get('species'):
                taxonomy.species = row['species']
        else:
            # Create new taxonomy
            taxonomy = Taxonomy(
                taxID=str(taxid),
                order=row.get('order'),
                family=row.get('family'),
                genus=row.get('genus'),
                species=row.get('species')
            )
            session.add(taxonomy)
            session.flush()
            logger.info(f"Created new taxonomy: {taxid}")
        
        return taxonomy
    
    def update_family_names(self, session, row: Dict[str, Any]) -> Optional[FamilyNames]:
        """Update or create family names record"""
        family_name = row.get('familyName')
        if not family_name:
            return None
            
        family = session.query(FamilyNames).filter(
            FamilyNames.familyName == family_name
        ).first()
        
        if family:
            # Update existing family
            if row.get('type_element_reference'):
                family.type_element_reference = row['type_element_reference']
        else:
            # Create new family
            family = FamilyNames(
                familyName=family_name,
                type_element_reference=row.get('type_element_reference')
            )
            session.add(family)
            session.flush()
            logger.info(f"Created new family: {family_name}")
        
        return family
    
    def get_navis_id(self, session, navis_name: str, family_id: Optional[int] = None) -> Optional[int]:
        """Get or create navis_names record and return its ID"""
        if not navis_name:
            return None
            
        # Query using ORM
        navis = session.query(NavisNames).filter(
            NavisNames.navis_name == navis_name
        ).first()
        
        if navis:
            return navis.id
        
        # Create new navis record
        navis = NavisNames(
            navis_name=navis_name,
            ship_family_id=family_id
        )
        session.add(navis)
        session.flush()  # Get the ID
        logger.info(f"Created new navis: {navis_name}")
        return navis.id
    
    def get_haplotype_id(self, session, haplotype_name: str, navis_id: Optional[int] = None, family_id: Optional[int] = None) -> Optional[int]:
        """Get or create haplotype_names record and return its ID"""
        if not haplotype_name:
            return None
            
        # Query using ORM
        haplotype = session.query(HaplotypeNames).filter(
            HaplotypeNames.haplotype_name == haplotype_name
        ).first()
        
        if haplotype:
            return haplotype.id
        
        # Create new haplotype record
        haplotype = HaplotypeNames(
            haplotype_name=haplotype_name,
            navis_id=navis_id,
            ship_family_id=family_id
        )
        session.add(haplotype)
        session.flush()  # Get the ID
        logger.info(f"Created new haplotype: {haplotype_name}")
        return haplotype.id
    
    def get_or_create_ship(self, session, accession: Accessions, row: Dict[str, Any]) -> Optional[Ships]:
        """Get or create a ship record for the given accession"""
        # Look for existing ship for this accession
        ship = session.query(Ships).filter(
            Ships.accession_id == accession.id
        ).first()
        
        if ship:
            return ship
        
        # Create new ship record - we don't have sequence data in TSV, so create minimal record
        ship = Ships(
            accession_id=accession.id,
            # We could add other fields here if available in the TSV
        )
        session.add(ship)
        session.flush()  # Get the ID
        logger.info(f"Created new ship for accession: {accession.accession_tag}")
        return ship
    
    def update_joined_ships(self, session, row: Dict[str, Any], accession: Accessions, 
                          family: Optional[FamilyNames]) -> bool:
        """Update or create joined ships record"""
        accession_tag = row.get('accession_tag')
        starship_id = row.get('starshipID')
        
        if not accession_tag or not starship_id:
            logger.warning(f"Missing required fields for joined_ships: {accession_tag}/{starship_id}")
            return False
        
        # Get or create the ship record for this accession
        ship = self.get_or_create_ship(session, accession, row)
        if not ship:
            logger.warning(f"Could not get/create ship for accession: {accession_tag}")
            return False
        
        # Get navis and haplotype IDs
        family_id = family.id if family else None
        navis_id = self.get_navis_id(session, row.get('starship_navis'), family_id)
        haplotype_id = self.get_haplotype_id(session, row.get('starship_haplotype'), navis_id, family_id)
        
        # Get taxonomy ID
        tax_id = None
        if row.get('taxID'):
            taxonomy = session.query(Taxonomy).filter(
                Taxonomy.taxID == str(row['taxID'])
            ).first()
            if taxonomy:
                tax_id = taxonomy.id
        
        # Find existing joined_ships record by starshipID 
        # Note: Using starshipID as primary identifier since that's what seems to be unique
        joined_ship = session.query(JoinedShips).filter(
            JoinedShips.starshipID == starship_id
        ).first()
        
        # Prepare the data
        ship_data = {
            'starshipID': starship_id,
            'curated_status': row.get('curated_status'),
            'ship_id': ship.id,
            'ship_family_id': family_id,
            'tax_id': tax_id,
            'ship_navis_id': navis_id,
            'ship_haplotype_id': haplotype_id,
        }
        
        if joined_ship:
            # Update existing record
            for key, value in ship_data.items():
                if hasattr(joined_ship, key):
                    setattr(joined_ship, key, value)
            logger.info(f"Updated joined_ships record: {accession_tag}/{starship_id}")
        else:
            # Create new record
            joined_ship = JoinedShips(**ship_data)
            session.add(joined_ship)
            logger.info(f"Created new joined_ships record: {accession_tag}/{starship_id}")
        
        return True
    
    def process_row(self, session, row: Dict[str, Any]) -> bool:
        """Process a single row from TSV"""
        try:
            self.stats['processed'] += 1
            
            # Update related tables first
            accession = self.update_accession(session, row)
            if not accession:
                logger.warning(f"No accession found/created for row {self.stats['processed']}")
                self.stats['skipped'] += 1
                return False
            
            taxonomy = self.update_taxonomy(session, row)
            family = self.update_family_names(session, row)
            
            # Update the main joined_ships table
            success = self.update_joined_ships(session, row, accession, family)
            
            if success:
                self.stats['updated'] += 1
                return True
            else:
                self.stats['skipped'] += 1
                return False
                
        except Exception as e:
            logger.error(f"Error processing row {self.stats['processed']}: {e}")
            self.stats['errors'] += 1
            return False
    
    def update_database(self, dry_run: bool = False) -> Dict[str, int]:
        """Main method to update database from TSV"""
        logger.info(f"Starting database update from {self.tsv_file_path}")
        logger.info(f"Dry run mode: {dry_run}")
        
        # Read TSV data
        tsv_data = self.read_tsv_data()
        
        with get_starbase_session() as session:
            try:
                for row in tsv_data:
                    try:
                        self.process_row(session, row)
                        
                        # Commit every 10 rows to avoid large transactions and handle errors per row
                        if self.stats['processed'] % 10 == 0:
                            if not dry_run:
                                session.commit()
                            logger.info(f"Processed {self.stats['processed']} rows...")
                    except Exception as row_error:
                        logger.error(f"Error processing row {self.stats['processed']}: {row_error}")
                        self.stats['errors'] += 1
                        session.rollback()  # Rollback failed transaction
                        continue  # Continue with next row
                
                # Final commit
                if not dry_run:
                    session.commit()
                    logger.info("Database update completed successfully")
                else:
                    session.rollback()
                    logger.info("Dry run completed - no changes made to database")
                
            except Exception as e:
                session.rollback()
                logger.error(f"Database update failed: {e}")
                raise
        
        return self.stats


def main():
    """Main function"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Update SQLite database from TSV file')
    parser.add_argument('tsv_file', help='Path to the TSV file')
    parser.add_argument('--dry-run', action='store_true', 
                       help='Perform a dry run without making changes')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Setup logging level
    if args.verbose:
        import logging
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Validate TSV file exists
    if not os.path.exists(args.tsv_file):
        print(f"Error: TSV file not found: {args.tsv_file}")
        sys.exit(1)
    
    try:
        # Create updater and run
        updater = DatabaseUpdater(args.tsv_file)
        stats = updater.update_database(dry_run=args.dry_run)
        
        # Print summary
        print("\n" + "="*50)
        print("DATABASE UPDATE SUMMARY")
        print("="*50)
        print(f"Processed: {stats['processed']}")
        print(f"Updated:   {stats['updated']}")
        print(f"Skipped:   {stats['skipped']}")
        print(f"Errors:    {stats['errors']}")
        print("="*50)
        
        if stats['errors'] > 0:
            print("⚠️  Some errors occurred. Check the logs above.")
            sys.exit(1)
        else:
            print("✅ Update completed successfully!")
            
    except Exception as e:
        print(f"❌ Fatal error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
