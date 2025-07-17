# Classification Workflow Testing

This directory contains scripts for testing and validating the classification workflow against known classifications in the database.

## Overview

The classification workflow testing system allows you to:
- Test the classification pipeline against all sequences in your database
- Compare predicted classifications with actual/known classifications
- Generate detailed performance metrics and confusion matrices
- Identify areas where the classification workflow needs improvement

## Files

- `test_classification_workflow.py` - Main testing framework
- `run_classification_tests.py` - Interactive script with example usage
- `TESTING_README.md` - This documentation file

## How It Works

The testing system works by:

1. **Fetching sequences** from the database along with their known classifications
2. **Running each sequence** through the classification workflow pipeline
3. **Comparing results** between predicted and actual classifications
4. **Generating metrics** including accuracy, confusion matrices, and detailed reports

## Classification Pipeline Stages

The workflow tests each sequence through these stages in order:

1. **Exact Match** - MD5 hash comparison for identical sequences
2. **Contained Match** - Check if sequence is contained within existing sequences
3. **Similar Match** - k-mer similarity comparison
4. **Family Classification** - HMMER-based family assignment
5. **Navis Classification** - Captain sequence clustering
6. **Haplotype Classification** - Sequence similarity within navis groups

## Usage

### Command Line Interface

```bash
# Run quick test (10 sequences)
python test_classification_workflow.py --max-sequences 10

# Test curated sequences only  
python test_classification_workflow.py --curated-only --max-sequences 50

# Full test with custom output directory
python test_classification_workflow.py --output-dir my_test_results

# Verbose output
python test_classification_workflow.py --verbose --max-sequences 5
```

#### Command Line Options

- `--output-dir` - Directory to save results (default: "classification_test_results")
- `--curated-only` - Only test curated sequences
- `--max-sequences` - Maximum number of sequences to test
- `--verbose` - Enable verbose logging

### Interactive Interface

```bash
python run_classification_tests.py
```

This provides a menu-driven interface with options for:
1. Quick test (10 sequences)
2. Curated sequences test (50 sequences)  
3. Full database test (all sequences)
4. Analyze previous results
5. Exit

### Programmatic Usage

```python
from test_classification_workflow import ClassificationTester

# Create tester instance
tester = ClassificationTester(
    output_dir="my_results",
    curated_only=True,
    max_sequences=100
)

# Run tests
tester.run_tests()
```

## Output Files

The testing system generates several output files:

### 1. Detailed Results (`detailed_results_TIMESTAMP.csv`)
Contains one row per tested sequence with columns:
- `accession_tag` - Sequence identifier
- `sequence_length` - Length of sequence
- `actual_family`, `actual_navis`, `actual_haplotype` - Known classifications
- `match_stage` - Which stage made the final classification
- `predicted_family`, `predicted_navis`, `predicted_haplotype` - Predicted classifications
- `found_match` - Whether any classification was found
- `family_correct`, `navis_correct`, `haplotype_correct` - Correctness for each stage
- `any_correct` - Whether any classification was correct
- `family_status`, `navis_status`, `haplotype_status` - Status of each classification stage:
  - `predicted` - Stage ran and made a prediction
  - `failed` - Stage ran but returned no result
  - `skipped` - Stage was skipped due to earlier stage failure
  - `error` - Stage encountered an error
  - `not_tested` - Stage was never attempted

### 2. Metrics Summary (`metrics_TIMESTAMP.json`)
JSON file containing:
- Overall accuracy metrics
- Classification counts by stage
- Accuracy by classification type
- Error counts

### 3. Confusion Matrices (`confusion_matrix_TYPE_TIMESTAMP.csv`)
Separate confusion matrices for:
- Family classification
- Navis classification  
- Haplotype classification

### 4. Summary Report (`summary_report_TIMESTAMP.txt`)
Human-readable summary including:
- Overall performance metrics
- Classification breakdown by stage
- Accuracy by classification type
- Confusion matrices

## Example Output

```
TEST SUMMARY
==================================================
Total sequences: 100
Classified sequences: 78
Non-sequences (exact/contained/similar): 22
Classification accuracy: 85.20%

Stage Status Summary:
  Family: 78 predicted, 12 failed, 0 skipped
  Navis: 45 predicted, 33 failed, 12 skipped  
  Haplotype: 23 predicted, 22 failed, 45 skipped

Results saved to: classification_test_results
```

## Interpreting Results

### Accuracy Metrics

- **Overall Accuracy** - Percentage of sequences where ANY classification (family, navis, haplotype, or exact match) was correct
- **Family Accuracy** - Percentage of family classifications that were correct
- **Navis Accuracy** - Percentage of navis classifications that were correct  
- **Haplotype Accuracy** - Percentage of haplotype classifications that were correct

### Classification Stages

Results show which stage of the pipeline made successful classifications:
- `exact` - Exact sequence matches (highest confidence)
- `contained` - Sequence contained in longer sequence
- `similar` - High similarity matches
- `family` - HMMER family classification
- `navis` - Captain clustering
- `haplotype` - Haplotype assignment

### Confusion Matrices

Show the relationship between actual vs predicted classifications:
- Diagonal values = correct predictions
- Off-diagonal values = misclassifications
- Can identify systematic errors or biases

## Performance Considerations

- **Quick tests** (10-50 sequences) - Complete in minutes
- **Curated tests** (100-500 sequences) - May take 30 minutes to 2 hours
- **Full database tests** (1000+ sequences) - Can take many hours or days

The time depends on:
- Database size
- Sequence lengths
- Available computational resources
- External tool performance (MMseqs2, BLAST, etc.)

## Troubleshooting

### Common Issues

1. **Database Connection Errors**
   - Ensure database is accessible and credentials are correct
   - Check that required tables exist

2. **Missing Dependencies**
   - Install required bioinformatics tools (MMseqs2, BLAST, etc.)
   - Ensure all Python dependencies are available

3. **Memory Issues**
   - Reduce `max_sequences` for large databases
   - Monitor system resources during testing

4. **Slow Performance**
   - Use `--curated-only` to focus on high-quality data
   - Limit sequences with `--max-sequences`
   - Run on systems with more CPU cores

### Debug Mode

Enable verbose logging to see detailed information:

```bash
python test_classification_workflow.py --verbose --max-sequences 5
```

This will show:
- Database queries being executed
- Classification stages being processed
- Detailed error messages
- Timing information

## Best Practices

1. **Start Small** - Begin with quick tests to verify the system works
2. **Test Curated Data** - Focus on high-quality, manually curated sequences first
3. **Monitor Resources** - Watch CPU and memory usage during large tests
4. **Save Results** - Keep detailed results for comparison after workflow improvements
5. **Iterate** - Use results to identify and fix classification issues

## Extending the Tests

The testing framework is designed to be extensible. You can:

1. **Add new metrics** by modifying `calculate_metrics()`
2. **Add new classification comparisons** in `compare_classifications()`
3. **Add new output formats** in the `save_results()` method
4. **Customize the workflow** by modifying `run_single_test()`

## Related Files

- `src/utils/classification_utils.py` - Main classification workflow
- `src/database/sql_manager.py` - Database access functions
- `src/config/logging.py` - Logging configuration 