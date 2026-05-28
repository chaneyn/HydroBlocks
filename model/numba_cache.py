import os
from pathlib import Path


def _iter_numba_cache_files(root_path):

 root_path = Path(root_path).expanduser()
 if not root_path.exists():
  return

 for suffix in ('*.nbc', '*.nbi', '*.pyc'):
  for cache_file in root_path.rglob(suffix):
   yield cache_file


def clear_numba_cache(comm, metadata, cache_roots, label):

 clear_at_start = metadata.get('numba_cache', {}).get('clear_at_start', False)
 if clear_at_start is False:
  return False

 roots = list(cache_roots)
 numba_cache_dir = os.environ.get('NUMBA_CACHE_DIR', '')
 if numba_cache_dir:
  roots.append(numba_cache_dir)

 rank = comm.Get_rank()
 if rank == 0:
  removed = 0
  for root in roots:
   for cache_file in _iter_numba_cache_files(root):
    try:
     cache_file.unlink()
     removed += 1
    except FileNotFoundError:
     pass
  print(f'{label}: cleared {removed} numba cache files before startup', flush=True)

 comm.Barrier()
 return True