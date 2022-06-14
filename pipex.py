import os
import sys
import datetime
import psutil
import fnmatch


def batch_processor():
    batch_filename = './pipex_batch_list.txt'
    python_command = './bin/python -u '
    pidfile_filename = './RUNNING'
    if "PIPEX_DATA" not in os.environ:
        os.environ['PIPEX_DATA'] = './data'
    if "PIPEX_WORK" in os.environ:
        batch_filename = './work/pipex_batch_list.txt'
        python_command = 'python -u '
        pidfile_filename = './work/RUNNING'

    print(">>> Start time pipex =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

    swap_used = False
    batch_file = open(batch_filename, 'r')
    while True:
        try:
            with open(pidfile_filename,'r') as f:
                lines = f.readlines()
                if psutil.pid_exists(int(lines[0])):
                    print(">>> Another PIPEX process seems to be running, exiting =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)
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
        elif curr_command.startswith('swap'):
            #creating required swap via bash script 'enable_swap.sh'
            print(">>> Creating swap space =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)
            swap_req = curr_command.replace('swap', '').strip()
            os.system('./enable_swap.sh ' + swap_req)
            swap_used = True
        elif len(curr_command.strip()) > 5:
            #pipex command
            print(">>> Processing next job =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)
            print('>>>    ' + curr_command.strip())
            os.system(python_command + curr_command.strip())
 
    batch_file.close()

    if "PIPEX_WORK" in os.environ:
        os.system('chmod -R 777 /opt/pipex/work')

    if swap_used:
        #deleting previously created swap via bash script 'disable_swap.sh'
        print(">>> Deleting swap space =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)
        os.system('./disable_swap.sh')

    for file in os.listdir('.'):
        if fnmatch.fnmatch(file, '<property*'):
            os.rmdir(file)
    
    if os.path.exists(pidfile_filename):
        os.remove(pidfile_filename)

    print(">>> End time pipex =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)


if __name__ == '__main__':
    batch_processor()
