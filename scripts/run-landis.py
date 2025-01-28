#libraries
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
import argparse
#batch_file = "historic-ncar.bat"

# function to run landis batch files in the docker container
def run_batch_file_in_container(batch_file, mount_path, container_name):
    try:
        #cmd /c bash or wine cmd /c
        command = f'sudo docker exec -it {container_name} cmd /c {mount_path}/{batch_file}'
        subprocess.run(command, shell=True, check=True)
        print(f"Completed {batch_file}.")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
        return False

# function to process a replicate and control output
def process_replicate(rep_id, batch_files, mount_path, output_folder, container_name):
    local_output_folder = os.path.join(output_folder, f"rep{rep_id}")
    os.makedirs(local_output_folder, exist_ok=True)

    for scenario_name, batch_file in batch_files:
        local_batch_file = os.path.join(local_output_folder, batch_file)
        if not os.path.exists(local_batch_file):
            print(f"Missing {batch_file} for replicate {rep_id}.")
            continue
        
        if not run_batch_file_in_container(batch_file, mount_path, container_name):
            print(f"Failed {batch_file} for replicate {rep_id}.")
            return

        print(f"Scenario {scenario_name}, replicate {rep_id} done.")

#main processing function
def main(num_replicates, batch_files, mount_path, output_folder, container_name):
    with ThreadPoolExecutor() as executor:
        executor.map(lambda rep_id: process_replicate(rep_id, batch_files, mount_path, output_folder, container_name), range(1, num_replicates + 1))

#using arg parser to take specified environment variables from the cli and use as vars in this script
if __name__ == '__main__':
    # command line arg parsing
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("--num_replicates", type=int, required=True, help="Number of replicates to process")
    arg_parser.add_argument("--batch_files", nargs='+', type=str, 
                            default=['historic-ncar.bat', 'future-ncar.bat', 'future-gfdl.bat'], 
                            help="Batch files to run (comma-separated)")
    arg_parser.add_argument("--container_name", type=str, required=True, help="Docker container name")
    arg_parser.add_argument("--access_id", type=str, required=True, help="User access ID")

    #parse
    args = arg_parser.parse_args()

    #split the batch files if many
    #batch_files = [(name.strip(), name.strip() + '.sh') for name in args.batch_files]

    # get environment vars
    num_replicates = args.num_replicates
    batch_files = args.batch_files
    print(f"Number of reps: {num_replicates}")
    print(f"Scenarios: {batch_files}")

    access_id = args.access_id
    container_name = args.container_name
    print(f"Access ID: {access_id}")
    print(f"Container name: {container_name}")

    # set mount path and output folder
    mount_path = f'/home/{access_id}/landis/input/' #mounting locally in js2 instance
    output_folder = f'/home/{access_id}/landis/output/'
    print(f"Output folder: {output_folder}")

    # run the main function
    main(num_replicates, mount_path, output_folder, container_name)
