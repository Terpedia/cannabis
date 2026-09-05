"""Unreviewed exact-family follow-up; NAD and NADP equations stay separate."""
from .phase1_missing_reference_review import run


if __name__ == '__main__':
    run(prefix='phase1-ureide-plant-reference',
        followup_statuses={'no-reference-sequence', 'weak-hits-only', 'no-hits'},
        source_prefix='phase1-ureide-gap')
