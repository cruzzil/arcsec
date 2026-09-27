# Security policy

## Supported versions

Security fixes are made to the latest release only.

## Reporting a vulnerability

Please report vulnerabilities privately through GitHub's
[private vulnerability reporting](https://github.com/cruzzil/arcsec/security/advisories/new)
rather than in a public issue. Include the arcsec version, the platform, and the steps
or input file that trigger the problem.

arcsec parses untrusted image files (FITS, XISF, ASDF) and catalogue files, and
downloads catalogues over HTTPS, so problems in those areas - crashes or excessive
memory use on a malformed file, path traversal when unpacking a catalogue archive -
are in scope.

You should receive a reply within a week. Fixes are released as a new patch version and
noted in `CHANGELOG.md`.
