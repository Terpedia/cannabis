# Cannabis service account

Created and verified on 2026-09-04.

- Project: `terpedia-489015`
- Principal: `cannabis-metabolome@terpedia-489015.iam.gserviceaccount.com`
- Project role: `roles/bigquery.jobUser`
- Dataset: `terpedia-489015.terpedia_core`
- Dataset role: `projects/terpedia-489015/roles/CannabisMetabolomeData`

The custom dataset role permits `bigquery.datasets.get` and
`bigquery.tables.create`, `get`, `getData`, `list`, `update`, and `updateData`.
It does not grant deletion or IAM administration. Write permissions can still
replace table contents, so ingestion must continue to use new, versioned tables
and check destination existence first.

The private key is stored outside this repository at
`/Users/danielmcshan/.config/gcloud/cannabis-metabolome/key.json`, with mode
`0600` inside a `0700` directory. Never commit, publish, or print the key.
The existing user login and Application Default Credentials were not replaced.

For local pipeline commands, set both variables in the invoking shell:

```sh
export GOOGLE_APPLICATION_CREDENTIALS=/Users/danielmcshan/.config/gcloud/cannabis-metabolome/key.json
export CLOUDSDK_AUTH_CREDENTIAL_FILE_OVERRIDE="$GOOGLE_APPLICATION_CREDENTIALS"
```

Use project `terpedia-489015` and BigQuery location `us-central1` explicitly.
Client libraries use the first variable; `gcloud` and `bq` use the second.
Cloud-hosted jobs should attach this service account as their runtime identity
instead of copying the local private key.

Verification: a BigQuery query returned `SESSION_USER()` equal to the service
account above and counted 268,924 rows in `terpene_identity_set`.
Subsequently, four pending versioned evidence tables were loaded under this
identity: archived references (161 rows), archived protein search (270 rows),
archived evidence (67 rows), and completion connectivity (7,396 rows).
Each destination was absent before loading. Job identities, destinations and
successful completion were verified, and every complete stored record matched
its local export. Receipts are `data/reports/phase1-archived-gcp.json` and
`data/reports/phase1-completion-connectivity-gcp.json`.

The dataset grant was made additively using SQL `GRANT`, preserving existing
access entries. The CLI dataset IAM-binding operation required allowlisting;
the SQL grant succeeded without broadening project-wide data access.
