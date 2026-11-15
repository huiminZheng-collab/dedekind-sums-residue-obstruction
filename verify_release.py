"""Verify this frozen release; mathematical replay uses a temporary copy."""
import argparse
import hashlib
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import sys
import tempfile

KEY = '490689A7D9DB559372C8620192664F521ACFE6E2'
ROOT_CERT = '552F7BDCF1A7AF9E6CE672017F4F12ABF77240C78E761AC203D1D9D20AC89988'
ENVELOPE = {'RELEASE-MANIFEST.sha256', 'RELEASE-MANIFEST.sha256.asc',
            'timestamp/manifest-signature.tsq', 'timestamp/manifest-signature.tsr',
            'timestamp/VERIFICATION.txt', 'artifacts/reproducibility-v1.0.0.zip'}
COMMANDS = [('src/verify.py',),('src/check_words.py','certificates/gamma1_9_words.txt'),
            ('src/verify_chi34.py',),('src/check_words.py','certificates/gamma1_12_words.txt'),
            ('src/audit_values.py',),('src/verify_k7.py',),
            ('src/verify_refinements.py',),('src/verify_k6.py',)]

def require(ok,message):
    if not ok:
        raise RuntimeError(message)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def run(args,cwd=None):
    p=subprocess.run([str(a) for a in args],cwd=cwd,text=True,
                     stdout=subprocess.PIPE,stderr=subprocess.STDOUT,encoding='utf8',errors='replace')
    require(p.returncode==0,'Command failed: '+str(args[0])+'\n'+p.stdout)
    return p.stdout

def gpath(path,gpg):
    s=str(Path(path).resolve())
    if os.name=='nt' and '/git/usr/bin/' in str(gpg).replace('\\','/').lower():
        return '/'+s[0].lower()+s[2:].replace('\\','/')
    return s

def check_hashes(root):
    expected={}
    for line in (root/'RELEASE-MANIFEST.sha256').read_text(encoding='ascii').splitlines():
        m=re.fullmatch(r'([a-f0-9]{64})  (.+)',line)
        require(m is not None,'Malformed manifest line')
        digest,name=m.groups();p=PurePosixPath(name)
        require(not p.is_absolute() and '..' not in p.parts and '\\' not in name,'Unsafe manifest path')
        require(name not in expected and name not in ENVELOPE,'Duplicate/circular manifest entry')
        expected[name]=digest
        path=root.joinpath(*p.parts)
        require(path.is_file() and not path.is_symlink(),'Missing/unsafe payload: '+name)
        require(sha(path)==digest,'Payload hash mismatch: '+name)
    actual={p.relative_to(root).as_posix() for p in root.rglob('*')
            if p.is_file() and '.git' not in p.relative_to(root).parts}
    require(not (actual-set(expected)-ENVELOPE),'Unexpected files: '+str(sorted(actual-set(expected)-ENVELOPE)))
    print('SHA-256 payload checks PASS:',len(expected),'files',flush=True)

def check_crypto(root,gpg,openssl):
    gpg=shutil.which(gpg) or gpg;openssl=shutil.which(openssl) or openssl
    with tempfile.TemporaryDirectory(prefix='ddk-') as folder:
        home=Path(folder)/'keys';home.mkdir()
        gp=lambda p:gpath(p,gpg)
        base=[gpg,'--homedir',gp(home),'--batch','--no-auto-key-retrieve']
        run(base+['--import',gp(root/'provenance/release-key.asc')])
        status=run(base+['--status-fd','1','--verify',gp(root/'RELEASE-MANIFEST.sha256.asc'),gp(root/'RELEASE-MANIFEST.sha256')])
        require('[GNUPG:] VALIDSIG '+KEY+' ' in status,'Unexpected signing-key fingerprint')
    cert=root/'provenance/DigiCertTrustedRootG4.pem'
    fp=run([openssl,'x509','-in',cert,'-noout','-fingerprint','-sha256']).split('=')[-1].strip().replace(':','')
    require(fp.upper()==ROOT_CERT,'Untrusted TSA root fingerprint')
    common=[openssl,'ts','-verify','-in',root/'timestamp/manifest-signature.tsr',
            '-CAfile',cert,'-untrusted',root/'provenance/DigiCertTrustedG4TimeStampingRSA4096SHA2562025CA1.pem']
    run(common+['-data',root/'RELEASE-MANIFEST.sha256.asc'])
    run(common+['-queryfile',root/'timestamp/manifest-signature.tsq'])
    print('OpenPGP signature, pinned TSA chain and both RFC 3161 imprint checks PASS',flush=True)

def replay(root):
    with tempfile.TemporaryDirectory(prefix='ddk-replay-') as folder:
        dest=Path(folder)/'release';shutil.copytree(root,dest,ignore=shutil.ignore_patterns('.git'))
        for args in COMMANDS:
            out=run([sys.executable,'-B',*args],cwd=dest)
            require('Traceback' not in out,'Python error in replay')
            print('Replay PASS:', ' '.join(args),flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--gpg',default='gpg');p.add_argument('--openssl',default='openssl');p.add_argument('--replay',action='store_true')
    a=p.parse_args();root=Path(__file__).resolve().parent
    check_hashes(root);check_crypto(root,a.gpg,a.openssl)
    if a.replay: replay(root)
    print('RELEASE VERIFICATION OK',flush=True)

if __name__=='__main__':
    main()
