import os
import sys
import psutil
from pipex_utils import log
import shutil
import fnmatch
import random
import subprocess


def batch_processor():
    batch_filename = './pipex_batch_list.txt'
    if os.path.exists("./Scripts/python.exe"):
        python_command = './Scripts/python.exe -u -W ignore '
    elif os.path.exists("./bin/python"):
        python_command = './bin/python -u -W ignore '
    else:
        python_command = 'python -u -W ignore '
    pidfile_filename = './RUNNING'
    run_id_filename = './run_id.txt'
    run_id_result_filename = './LAST_RUN_ID'
    run_id = random.randint(1, 1000000)
    log_filename = './log.txt'
    if "PIPEX_DATA" not in os.environ:
        os.environ['PIPEX_DATA'] = './data'
    if "PIPEX_WORK" in os.environ:
        batch_filename = './work/pipex_batch_list.txt'
        python_command = 'python -u '
        pidfile_filename = './work/RUNNING'
        run_id_filename = './work/run_id.txt'
        run_id_result_filename = './work/LAST_RUN_ID'
        log_filename = './work/log.txt'

    log("Start time pipex")

    if os.path.exists(run_id_filename):
        with open(run_id_filename, 'r', encoding='utf-8') as f:
            run_id = f.read().strip()
            log("Setting run_id = " + str(run_id))

    swap_used = False
    batch_file = open(batch_filename, 'r')

    try:
        while True:
            try:
                with open(pidfile_filename,'r') as f:
                    lines = f.readlines()
                    if psutil.pid_exists(int(lines[0])):
                        log("Another PIPEX process seems to be running, exiting")
                        sys.exit()
            except IOError:
                pass

            curr_command = batch_file.readline()
            if not curr_command:
                #EOF
                break
            elif curr_command.startswith('#'):
                #comment
                continue
            elif curr_command.startswith('run_id'):
                #using the run_id provided by the user
                run_id = curr_command.replace('run_id', '').strip()
                log("Setting run_id = " + str(run_id))
            elif curr_command.startswith('max_res'):
                #setting PIPEX_MAX_RESOLUTION environment variable
                os.environ["PIPEX_MAX_RESOLUTION"] = curr_command.replace('max_res', '').strip()
                log("Setting max resolution to = " + os.environ["PIPEX_MAX_RESOLUTION"])
            elif curr_command.startswith('swap'):
                #creating required swap via bash script 'enable_swap.sh'
                log("Creating swap space")
                swap_req = curr_command.replace('swap', '').strip()
                os.system('./enable_swap.sh ' + swap_req)
                swap_used = True
            elif len(curr_command.strip()) > 5:
                #pipex command
                try:
                    log("Processing next job")
                    print('>>>    ' + curr_command.strip())
                    os.system(python_command + curr_command.strip())
                except Exception:
                    pass

                if '-data=' in curr_command:
                    arg_start_index = curr_command.index('-data=') + 6
                    end_char = ' '
                    if curr_command[arg_start_index:arg_start_index + 1] == '\'':
                        end_char = '\''
                        arg_start_index = arg_start_index + 1
                    elif curr_command[arg_start_index:arg_start_index + 1] == '\"':
                        end_char = '\"'
                        arg_start_index = arg_start_index + 1
                    arg_end_index = curr_command.index(end_char, arg_start_index + 1)
                    curr_data_folder = curr_command[arg_start_index:arg_end_index].strip()
                    if os.path.exists(log_filename):
                        shutil.copyfile(log_filename, os.path.join(curr_data_folder, os.path.basename(log_filename)))

        batch_file.close()

    except Exception as e:
        print(e)

    if "PIPEX_WORK" in os.environ:
        os.system('chmod -R 777 /opt/pipex/work')

    if swap_used:
        #deleting previously created swap via bash script 'disable_swap.sh'
        log("Deleting swap space")
        os.system('./disable_swap.sh')

    for file in os.listdir('.'):
        if fnmatch.fnmatch(file, '<property*'):
            os.rmdir(file)
    
    if os.path.exists(pidfile_filename):
        os.remove(pidfile_filename)

    with open(run_id_result_filename, 'w', encoding='utf-8') as f:
        f.write(str(run_id))
        f.close()

    log("End time pipex")


if __name__ == '__main__':
    batch_processor()
