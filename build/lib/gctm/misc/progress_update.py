import os
import sys

# Functions used to print progress of various processes


def progress_func_time_dep(func):

    return print("Calculated "+str(func))


def progress_func_power(i, nk):

    if nk <= 1000:

        if i == 0:

            print("Calculating the power spectrum:")

        else:

            n_bar = 10  # size of progress bar

            j = i/nk

            sys.stdout.write('\r')
            sys.stdout.write(f"[{'=' * int(n_bar * j):{n_bar}s}] {int(100 * j)}%")
            sys.stdout.flush()

        if i == nk-1:

            sys.stdout.write('\n')

    if nk > 1000:

        if i == 0:

            print("Calculating the power spectrum:")

        else:

            n_bar = 10  # size of progress bar

            j = i/nk

            sys.stdout.write('\r')
            sys.stdout.write(f"[{'=' * int(n_bar * j):{n_bar}s}] {int(100 * j)}%")
            sys.stdout.flush()


        if i == nk-1:

            sys.stdout.write('\n')    

    return


def progress_func_power_save(filename, i, nk):

    if i == 0:

        print("Started calculating power")

        with open(os.path.join("power_progress_"+filename+".txt"), 'a') as file:

            file.writelines("Started calculating power"+'\n')

    if i % 100 == 0:

        percent = (float(i)/float(nk))*100.0

        print("Completed "+str(percent)+"%")

        with open(os.path.join("power_progress_"+filename+".txt"), 'a') as file:

            file.writelines("Percent complete: "+str(percent)+"%"+'\n')

    return
