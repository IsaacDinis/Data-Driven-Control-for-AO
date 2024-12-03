import subprocess
import signal
import os
import time
from scipy.io import loadmat
import dao
import numpy as np
import astropy.io.fits as pfits
# List of script commands
scripts = [
    "daoPapyrusCreateShm.py",
    "daoPapyrusCreateSHM.py",
    "daoPapyrusSimStart",
    "papyCtrl.py",
    "daoDmDisp.py -s dmCmd -m dm241Map",
    "daoBarRTD.py -s /tmp/papyrus_modes.im.shm",
]

# List to store process objects
processes = []

# Function to clean up all processes and tmux sessions
def cleanup():
    print("\nCleaning up...")
    # Terminate all the processes
    for proc in processes:
        try:
            print(f"Terminating {proc.pid}...")
            os.kill(proc.pid, signal.SIGTERM)  # Send SIGTERM to the process
            print(f"Process {proc.pid} terminated.")
        except Exception as e:
            print(f"Error terminating process {proc.pid}: {e}")

    # Kill all tmux sessions
    try:
        print("Killing all tmux sessions...")
        subprocess.run(["tmux", "kill-server"], check=True)
        print("All tmux sessions terminated.")
    except FileNotFoundError:
        print("tmux not found. Ensure it is installed and in the PATH.")
    except subprocess.CalledProcessError as e:
        print(f"Error while killing tmux sessions: {e}")

    print("Cleanup complete.")

# Run the scripts in the background
for script in scripts:
    try:
        print(f"Running {script}...")
        # Split the command into a list for Popen compatibility
        command = script.split()
        proc = subprocess.Popen(command)  # Start the script in the background
        processes.append(proc)  # Keep track of the process
        print(f"{script} started with PID {proc.pid}.")
    except FileNotFoundError:
        print(f"{script} not found. Ensure it is in the PATH or the current directory.")
time.sleep(1) # wait for SHM to be created
amp=10
turb = pfits.getdata("tilt_traj.fits")
dmTurb=dao.shm('/tmp/dmCmd03.im.shm')
M2V = dao.shm("/tmp/m2c.im.shm").get_data()
# Infinite loop to wait for Ctrl+C
try:
    print("Press Ctrl+C to terminate all jobs and tmux sessions...")
    while True:  
        time.sleep(1)
        # for k in range(turb.shape[0]):
        #     dmTurb.set_data(M2V[:,0]*turb[k].astype(np.float32)*amp)
        #     time.sleep(0.01)
        # for k in np.linspace(turb.shape[0]-1,0,turb.shape[0]):
        #     dmTurb.set_data(M2V[:,0]*turb[int(k)].astype(np.float32)*amp)
        #     time.sleep(0.01)


except KeyboardInterrupt:
    # Handle Ctrl+C gracefully
    print("\nCtrl+C detected!")
    cleanup()
