# Joined Ships Duplicate Removal

## Overview

This utility removes duplicate entries in the `joined_ships` table that have `_ACC_XXXX` or `_SHIP_XXXX` suffixes in their `starshipID` column.

## What it does

1. **Identifies duplicates**: Finds entries with suffixes like `_ACC_1468` or `_SHIP_453` appended to their `starshipID`
2. **Retains originals**: Keeps the entry with the lowest ID (usually the first created)
3. **Coalesces data**: Before deleting duplicates, it copies any non-NULL values from the duplicate entries to fill NULL fields in the original entry
4. **Deletes duplicates**: Removes the duplicate entries

## Fields that are coalesced

The following fields will have their NULL values filled from duplicate entries:
- `accession_id`
- `ship_id`
- `ship_family_id`
- `tax_id`
- `genome_id`
- `captain_id`
- `ship_navis_id`
- `ship_haplotype_id`

## Usage

### Dry Run (Recommended First)

Always run a dry-run first to see what changes will be made:

```bash
conda activate starbase_minimal
python remove_joined_ships_duplicates.py --dry-run
```

This will:
- Show how many duplicate groups were found
- Display how many entries will be updated/deleted
- Show example coalesced values
- **Make NO changes to the database**

### Apply Changes

Once you've reviewed the dry-run output and are satisfied:

```bash
conda activate starbase_minimal
python remove_joined_ships_duplicates.py
```

### Save Detailed Report

To save a detailed JSON report of all changes:

```bash
conda activate starbase_minimal
python remove_joined_ships_duplicates.py --dry-run --output report.json
```

Or when applying changes:

```bash
conda activate starbase_minimal
python remove_joined_ships_duplicates.py --output report.json
```

## Example Output

```
================================================================================
SUMMARY
================================================================================
Duplicate groups found: 2826
Entries updated (coalesced): 1162
Entries deleted: 4423
Total fields coalesced: 1596

Example duplicate groups:

Group 1:
  Base starshipID: altals1_s00058
  Base ID: 1
  Duplicate count: 2
  Duplicate IDs: [1742, 3679]
  Fields coalesced: accession_id, ship_id

Update 1:
  Entry ID: 1
  StarshipID: altals1_s00058
  Fields updated:
    accession_id: NULL -> 453 (from altals1_s00058_SHIP_453)
    ship_id: NULL -> 453 (from altals1_s00058_SHIP_453)
```

## Database Schema

The script works with the `joined_ships` table which has the following structure:

```sql
CREATE TABLE joined_ships (
    id INTEGER PRIMARY KEY,
    starshipID TEXT,
    evidence TEXT,
    source TEXT,
    curated_status TEXT,
    accession_id INTEGER,
    ship_id INTEGER,
    ship_family_id INTEGER,
    tax_id INTEGER,
    genome_id INTEGER,
    captain_id INTEGER,
    ship_navis_id INTEGER,
    ship_haplotype_id INTEGER,
    created_at TIMESTAMP,
    updated_at TIMESTAMP
);
```

## Technical Details

### How Duplicates are Identified

The script uses a regular expression pattern to identify duplicates:
```python
r'^(.+?)_(ACC|SHIP)_(\d+)$'
```

This matches starshipIDs like:
- `08-19-02-30_s00028_ACC_1468` → Base: `08-19-02-30_s00028`, Type: `ACC`, Number: `1468`
- `altals1_s00058_SHIP_453` → Base: `altals1_s00058`, Type: `SHIP`, Number: `453`

### How Originals are Determined

For each group of duplicates:
1. If a base entry exists (without suffix), it's used as the original
2. Among multiple candidates, the one with the lowest `id` is chosen
3. This ensures the first-created entry is preserved

### Coalescing Strategy

For each NULL field in the original:
1. Check duplicate entries in order
2. Use the first non-NULL value found
3. Update the original entry with this value
4. Record the source of the value for tracking

## Safety Features

1. **Dry-run mode**: Test before applying changes
2. **Transaction rollback**: On any error, all changes are rolled back
3. **Detailed logging**: All actions are logged for audit trail
4. **JSON reports**: Optional detailed reports for review

## Troubleshooting

### ImportError: No module named 'sqlalchemy'

Make sure you've activated the correct conda environment:
```bash
conda activate starbase_minimal
```

### Database is locked

Make sure no other processes are accessing the database:
```bash
# Check for running processes
ps aux | grep python

# Stop the web application if running
# Then retry the script
```

### Script runs but no duplicates found

This means your database doesn't have entries with `_ACC_XXXX` or `_SHIP_XXXX` suffixes, which is good!

## Integration with Cleanup Pipeline

This functionality is also available as a function in the database cleanup utilities:

```python
from src.database.cleanup.utils.database_cleanup import remove_suffixed_joined_ships_duplicates

# Dry run
report = remove_suffixed_joined_ships_duplicates(dry_run=True)

# Apply changes
report = remove_suffixed_joined_ships_duplicates(dry_run=False)
```

## Related Scripts

- `src/database/cleanup/utils/database_cleanup.py` - Contains the core cleanup functions
- `src/database/cleanup/run_comprehensive_cleanup.py` - Comprehensive database cleanup pipeline

## Questions or Issues?

Check the logs in the console output for detailed information about what the script is doing.

