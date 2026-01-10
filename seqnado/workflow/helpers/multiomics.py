"""Helper functions for multiomics workflows."""


def get_assay_all_inputs(ASSAYS, rules):
    """
    Get all inputs from assay-specific 'all' rules.

    Args:
        ASSAYS: List of assay objects.
        rules: Snakemake rules object.

    Returns:
        list: List of inputs from all assay 'all' rules.
    """
    inputs = []
    for assay in ASSAYS:
        rule_name = f"{assay.clean_name}_all"
        inputs.append(getattr(rules, rule_name).input)
    return inputs
