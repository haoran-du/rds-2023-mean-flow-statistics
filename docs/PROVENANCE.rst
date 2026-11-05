Extraction provenance
=====================

This independent 2023 repository was extracted from the validated combined
workbench at commit 2bd418759368bdb305697ea354e62e99f257242f.
``extraction_provenance.json`` records source hashes and checks that the
computational function/class bodies are unchanged. Package imports, the
command-line interface and the paper-specific documentation were separated.

Small helpers are deliberately present locally, not imported from another
paper package. Each repository can be installed alone. Package namespaces
rds2023, rds2024 and rds2025 also permit installation together without a
shared ``rds`` namespace collision. No fourth shared-code repository is
required. Future shared numerical changes should be tested in each consumer.

The original workbench and its complete archive remain preserved. New
idealised-campaign work stays there rather than being represented as part
of this published paper. Experimental data have not been copied to GitHub.
