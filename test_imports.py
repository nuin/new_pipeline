import os
import sys
from pathlib import Path
import traceback

# Ensure we're in the project root
project_root = Path(__file__).parent
os.chdir(project_root)

print("Current working directory:")
print(Path.cwd())

# Add the project root to Python path
sys.path.insert(0, str(project_root))

print("\nPython path:")
for path in sys.path:
    print(path)

print("\nTrying to import src:")
try:
    import src
    print("Successfully imported src")
except ImportError as e:
    print(f"Failed to import src: {e}")
    print("Traceback:")
    traceback.print_exc()

print("\nTrying to import src.pipeline:")
try:
    from src import pipeline
    print("Successfully imported src.pipeline")
except ImportError as e:
    print(f"Failed to import src.pipeline: {e}")
    print("Traceback:")
    traceback.print_exc()

print("\nTrying to import src.lib:")
try:
    from src import lib
    print("Successfully imported src.lib")
except ImportError as e:
    print(f"Failed to import src.lib: {e}")
    print("Traceback:")
    traceback.print_exc()

print("\nTrying to import src.lib.bwa_align:")
try:
    from src.lib import bwa_align
    print("Successfully imported src.lib.bwa_align")
except ImportError as e:
    print(f"Failed to import src.lib.bwa_align: {e}")
    print("Traceback:")
    traceback.print_exc()

print("\nListing contents of src:")
try:
    print(dir(src))
except:
    print("Failed to list contents of src")

print("\nListing contents of src.lib:")
try:
    print(dir(src.lib))
except:
    print("Failed to list contents of src.lib")

print("\nListing contents of current directory:")
for item in Path.cwd().iterdir():
    print(item)

print("\nListing contents of src directory:")
src_dir = Path.cwd() / 'src'
if src_dir.exists():
    for item in src_dir.iterdir():
        print(item)
else:
    print("src directory not found")

print("\nListing contents of src/lib directory:")
lib_dir = src_dir / 'lib'
if lib_dir.exists():
    for item in lib_dir.iterdir():
        print(item)
else:
    print("src/lib directory not found")

print("\nChecking __init__.py files:")
init_files = [
    Path.cwd() / 'src' / '__init__.py',
    Path.cwd() / 'src' / 'lib' / '__init__.py'
]
for init_file in init_files:
    if init_file.exists():
        print(f"{init_file} exists")
        print("Contents:")
        print(init_file.read_text())
    else:
        print(f"{init_file} does not exist")