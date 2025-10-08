#!/usr/bin/env python3
"""
Script to link joined_ships to ships via FASTA sequence matching.

This script:
1. Parses a FASTA file containing starship sequences
2. Matches joined_ships.starshipID to FASTA headers (partial match)
3. Calculates MD5 hash of the sequence
4. Links to existing ship with matching MD5, or creates new ship
5. Updates joined_ships.ship_id with the link

This is the proper way to populate joined_ships.ship_id - based on actual
sequence identity rather than accession relationships.

Usage:
    python link_ships_via_fasta.py --fasta sequences.fasta --dry-run
    python link_ships_via_fasta.py --fasta sequences.fasta
    python link_ships_via_fasta.py --fasta sequences.fasta --output report.json
"""

import sys
import os
import argparse
import json
from datetime import datetime

# Add project root to path
PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from src.database.cleanup.utils.database_cleanup import link_joined_ships_via_fasta
from src.config.logging import get_logger

logger = get_logger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description='Link joined_ships to ships via FASTA sequence matching'
    )
    parser.add_argument(
        '--fasta',
        type=str,
        required=True,
        help='Path to FASTA file containing starship sequences'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Analyze and report what would be done without making changes'
    )
    parser.add_argument(
        '--output',
        type=str,
        help='Path to save the detailed report (JSON format)'
    )
    
    args = parser.parse_args()
    
    # Check if FASTA file exists
    if not os.path.exists(args.fasta):
        logger.error(f"FASTA file not found: {args.fasta}")
        return 1
    
    logger.info("=" * 80)
    logger.info("LINK JOINED_SHIPS VIA FASTA")
    logger.info("=" * 80)
    logger.info(f"Mode: {'DRY RUN' if args.dry_run else 'APPLYING CHANGES'}")
    logger.info(f"FASTA file: {args.fasta}")
    logger.info("")
    
    try:
        # Run the linking
        report = link_joined_ships_via_fasta(args.fasta, dry_run=args.dry_run)
        
        # Display summary
        logger.info("")
        logger.info("=" * 80)
        logger.info("SUMMARY")
        logger.info("=" * 80)
        logger.info(f"FASTA entries parsed: {report['summary']['fasta_entries']}")
        logger.info(f"Joined_ships processed: {report['summary']['joined_ships_processed']}")
        logger.info(f"Ships linked: {report['summary']['ships_linked']}")
        logger.info(f"Ships created: {report['summary']['ships_created']}")
        logger.info(f"Fallback FNA matches: {report['summary']['fallback_fna_matches']}")
        logger.info(f"No FASTA match: {report['summary']['no_fasta_match']}")
        logger.info(f"Warnings: {report['summary']['warnings']}")
        logger.info("")
        logger.info(f"Recommendation: {report['summary']['recommendation']}")
        logger.info("")
        
        # Show examples of ships linked
        if report['ships_linked']:
            logger.info("=" * 80)
            logger.info("Example ships linked (first 10):")
            for i, link in enumerate(report['ships_linked'][:10]):
                logger.info(f"{i+1}. JS {link['joined_ships_id']} ({link['starshipID']})")
                logger.info(f"   → Ship {link['ship_id']}, MD5: {link['md5'][:16]}...")
                logger.info(f"   FASTA: {link['fasta_header'][:60]}...")
        
        # Show examples of ships created
        if report['ships_created']:
            logger.info("")
            logger.info("=" * 80)
            logger.info(f"Ships created (total: {len(report['ships_created'])}):")
            logger.info("Examples (first 5):")
            for i, ship in enumerate(report['ships_created'][:5]):
                logger.info(f"{i+1}. Ship {ship['ship_id']} for {ship['starshipID']}")
                logger.info(f"   Sequence length: {ship['sequence_length']}, MD5: {ship['md5'][:16]}...")
                logger.info(f"   FASTA: {ship['fasta_header'][:60]}...")
        
        # Show fallback FNA matches
        if report['fallback_fna_matches']:
            logger.info("")
            logger.info("=" * 80)
            logger.info(f"Fallback FNA directory matches (total: {len(report['fallback_fna_matches'])}):")
            logger.info("Examples (first 10):")
            for i, match in enumerate(report['fallback_fna_matches'][:10]):
                logger.info(f"{i+1}. JS {match['joined_ships_id']} ({match['starshipID']})")
                logger.info(f"   Found in: {match['filename']}")
        
        # Show entries without FASTA match
        if report['no_fasta_match']:
            logger.info("")
            logger.info("=" * 80)
            logger.info(f"Entries without FASTA match: {len(report['no_fasta_match'])}")
            logger.info("Examples (first 10):")
            for i, nomatch in enumerate(report['no_fasta_match'][:10]):
                logger.info(f"{i+1}. JS {nomatch['joined_ships_id']}: {nomatch['starshipID']}")
                logger.info(f"   Issue: {nomatch['issue']}")
        
        # Show warnings
        if report['warnings']:
            logger.info("")
            logger.info("=" * 80)
            logger.info(f"Warnings: {len(report['warnings'])}")
            for i, warning in enumerate(report['warnings'][:10]):
                logger.info(f"{i+1}. JS {warning['joined_ships_id']} ({warning['starshipID']})")
                logger.info(f"   Issue: {warning['issue']}")
        
        # Save report if requested
        if args.output:
            output_path = args.output
            with open(output_path, 'w') as f:
                json.dump(report, f, indent=2, default=str)
            logger.info(f"\nDetailed report saved to: {output_path}")
        
        # Final message
        logger.info("")
        logger.info("=" * 80)
        if args.dry_run:
            logger.info("DRY RUN COMPLETE - No changes were made")
            logger.info("Run without --dry-run to apply these changes")
        else:
            logger.info("CHANGES APPLIED SUCCESSFULLY")
            logger.info("")
            logger.info("Summary:")
            logger.info(f"  - {report['summary']['ships_linked']} joined_ships entries linked to ships")
            logger.info(f"  - {report['summary']['ships_created']} new ships created from FASTA sequences")
            if report['summary']['fallback_fna_matches'] > 0:
                logger.info(f"  - {report['summary']['fallback_fna_matches']} matches found via fallback FNA directory")
            if report['summary']['no_fasta_match'] > 0:
                logger.info(f"  - {report['summary']['no_fasta_match']} entries had no FASTA match (need manual review)")
        logger.info("=" * 80)
        
        return 0
        
    except FileNotFoundError as e:
        logger.error(f"File not found: {str(e)}")
        return 1
    except Exception as e:
        logger.error(f"Error during FASTA linking: {str(e)}", exc_info=True)
        return 1


if __name__ == '__main__':
    sys.exit(main())

