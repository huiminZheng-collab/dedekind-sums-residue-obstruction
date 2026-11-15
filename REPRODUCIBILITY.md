# Reproducibility

The mathematical certificates require Python 3.12 and its standard library. No author repository, external data download, or symbolic-algebra package is required. Programs recompute finite sums with rational arithmetic. Their output files are written under `certificates/`; work on an extracted copy if the frozen release is to remain unchanged.

From the directory containing `src/` and `certificates/`, run in this order:

```text
python -B src/verify.py
python -B src/check_words.py certificates/gamma1_9_words.txt
python -B src/verify_chi34.py
python -B src/check_words.py certificates/gamma1_12_words.txt
python -B src/audit_values.py
python -B src/verify_k7.py
python -B src/verify_refinements.py
python -B src/verify_k6.py
```

The final script messages are respectively `CERTIFICATE OK`, `words ok:`, `CHI34 CERTIFICATE OK`, `words ok:`, `AUDIT VALUES OK`, `K7 CERTIFICATE OK`, `REFINEMENTS CERTIFICATE OK`, and `K6 CERTIFICATE OK`. The refinements program reads the weight-7 log, so the stated order matters.

These finite checks supply the input to the manuscript's cocycle and integrality proofs. They are not an enumeration of the infinite groups. Finite image scans mentioned as experiments in the paper remain experiments.

## Release provenance

To check hashes, the detached OpenPGP signature and RFC 3161 evidence, install GnuPG and OpenSSL and run:

```text
python -B verify_release.py --gpg gpg --openssl openssl
python -B verify_release.py --gpg gpg --openssl openssl --replay
```

Executable paths may be supplied instead of `gpg` or `openssl`. The verifier imports only the public release key into a temporary keyring. It does not request any private key or password, alter an existing keyring, or contact a network server. On Git for Windows, executable paths containing `Git/usr/bin` are handled as MSYS paths.

The release key has fingerprint `490689A7D9DB559372C8620192664F521ACFE6E2`. The pinned self-signed DigiCert Trusted Root G4 SHA-256 certificate fingerprint is `552F7BDCF1A7AF9E6CE672017F4F12ABF77240C78E761AC203D1D9D20AC89988`. The certificate was downloaded from DigiCert's official root-certificate directory and checked against the local trusted CA bundle before packaging.

The payload manifest covers the manuscript, source, code, saved finite certificates, public key, CA certificates and these instructions. It excludes its own signature and the later timestamp/archive envelope to avoid circular hashes. Those envelope objects are cryptographically verified separately. An unexpected non-envelope file is rejected; a local `.git/` directory is treated as version-control metadata.

The ZIP is constructed twice with fixed member order, timestamps and permissions and checked for identical bytes. This is a claim about packaging, not about byte-identical regeneration of the PDF by every TeX installation. The PDF is the frozen manuscript supplied with this release. To compile the source, use `pdflatex` twice in `paper/`; regenerated PDF metadata may differ.
