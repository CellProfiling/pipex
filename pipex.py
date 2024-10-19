import os
import sys
import datetime
import psutil
import shutil
import fnmatch
import random


def batch_processor():
    batch_filename = './pipex_batch_list.txt'
    python_command = 'python -u -W ignore '
    if os.path.exists("./bin/python"):
        python_command = './bin/python -u '
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

    print(">>> Start time pipex =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)

    if os.path.exists(run_id_filename):
        with open(run_id_filename, 'r', encoding='utf-8') as f:
            run_id = f.read().strip()
            print(">>> Setting run_id =",run_id, datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)

    swap_used = False
    batch_file = open(batch_filename, 'r')
    while True:
        try:
            with open(pidfile_filename,'r') as f:
                lines = f.readlines()
                if psutil.pid_exists(int(lines[0])):
                    print(">>> Another PIPEX process seems to be running, exiting =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
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
            print(">>> Setting run_id =",run_id, datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
        elif curr_command.startswith('swap'):
            #creating required swap via bash script 'enable_swap.sh'
            print(">>> Creating swap space =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
            swap_req = curr_command.replace('swap', '').strip()
            os.system('./enable_swap.sh ' + swap_req)
            swap_used = True
        elif len(curr_command.strip()) > 5:
            #pipex command
            try:
                print(">>> Processing next job =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
                print('>>>    ' + curr_command.strip())
                os.system(python_command + curr_command.strip())
            except Exception:
                pass

            if curr_command.index('-data=') > 0:
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
                    shutil.copyfile(log_filename, curr_data_folder + '/' + os.path.basename(log_filename))
 
    batch_file.close()

    if "PIPEX_WORK" in os.environ:
        os.system('chmod -R 777 /opt/pipex/work')

    if swap_used:
        #deleting previously created swap via bash script 'disable_swap.sh'
        print(">>> Deleting swap space =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
        os.system('./disable_swap.sh')

    for file in os.listdir('.'):
        if fnmatch.fnmatch(file, '<property*'):
            os.rmdir(file)
    
    if os.path.exists(pidfile_filename):
        os.remove(pidfile_filename)

    with open(run_id_result_filename, 'w', encoding='utf-8') as f:
        f.write(str(run_id))
        f.close()

    print(">>> End time pipex =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)


if __name__ == '__main__':
    batch_processor()
