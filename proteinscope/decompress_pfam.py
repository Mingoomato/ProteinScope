import gzip, shutil
from pathlib import Path

src = Path(r'C:\Users\jmkjm\Desktop\Protein_Finder_Project\proteinscope\.data\pfam\Pfam-A.hmm.gz')
dst = src.with_suffix('')

print('Decompressing', src.stat().st_size // 1_000_000, 'MB ...')
with gzip.open(src, 'rb') as f_in, open(dst, 'wb') as f_out:
    shutil.copyfileobj(f_in, f_out)
print('Done:', dst.stat().st_size // 1_000_000, 'MB uncompressed')
