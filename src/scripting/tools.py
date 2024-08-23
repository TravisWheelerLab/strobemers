import os

def experiment_dir(experiment_name):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    relative_path = f"../../data/{experiment_name}"
    relative_path_to_script = os.path.join(script_dir, relative_path)
    return relative_path_to_script