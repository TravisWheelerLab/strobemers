import sys
import os

"""
I want this script to automate experiments.
But I'm using a bash script instead :)
"""
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <directory_path>")
        sys.exit(1) 
    
    directory_path = sys.argv[1]
    assert os.path.isdir(directory_path)

