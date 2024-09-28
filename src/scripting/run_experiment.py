import sys
import os

""" I want this script to automate experiments.
1. Take a folder name.
2. In the folder, match query files to reference files.
3. Run strobemer program on pairs of query-reference files.
4. Combine outputted csv files from rust program into one.
"""
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <directory_path>")
        sys.exit(1)
    
    directory_path = sys.argv[1]
    assert os.path.isdir(directory_path)

