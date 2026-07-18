# Security Policy

## Supported versions

| Version | Supported |
| :-- | :--: |
| 0.1.x   | Yes |
| 0.0.x   | No  |

Only the latest minor version receives security and bug fixes.

## Reporting a vulnerability

If you discover a security issue in `geoeq`, please
[open a GitHub issue](https://github.com/geoeq/geoeq/issues/new) with the
label `security`.

We aim to:

1. Acknowledge the report within **72 hours**.
2. Triage and confirm within **7 days**.
3. Issue a patched release within **30 days** of confirmation, sooner for
   critical issues.

Please include in your report:

- A description of the vulnerability and the affected component.
- A minimal reproducer (Python script or test case).
- The version of `geoeq` and Python you tested with.
- Any suggested mitigation.

## Scope

`geoeq` is a pure-Python scientific library. Its security surface is small,
but the following categories of issue are in scope:

- **Code execution from untrusted input** — for example, a malicious AGS,
  GEF, or CSV file passed to `ge.read_ags`, `ge.read_gef`, or `ge.read_csv`
  that triggers arbitrary code execution, a crash that exposes data, or
  unsafe deserialisation.
- **Path traversal** in file readers.
- **Resource exhaustion** — pathological inputs that cause unbounded
  memory or CPU usage.
- **Dependency vulnerabilities** — security issues in `numpy`,
  `matplotlib`, or `scipy` that we should pin against.

The following are **out of scope**:

- Numerical inaccuracy or convergence issues in scientific formulas —
  please open a regular issue with a textbook citation.
- Compatibility issues with old Python versions or unsupported platforms.
- Issues in third-party tools that consume `geoeq` output.

## Responsible disclosure

We will credit security researchers in the release notes of the patched
version, unless requested otherwise. Please give us a reasonable window
to publish a fix before disclosing publicly.
