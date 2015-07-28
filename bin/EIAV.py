from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
# from timeit import default_timer
# from mpl_toolkits.mplot3d import Axes3D
import gc as garbage
import os
import json
import argparse
# import shlex
# from time import sleep
# from subprocess import Popen, PIPE
from datetime import datetime


DATE_TIME_DIRECTORY_FORMAT = '%y%m%d_%H%M%S'
DIRECTORY = os.path.dirname(os.path.abspath(__file__))
GRAPH_SAVED = """
GRAPH SAVED
-----------
%s
"""


def PARSE_ARGS():
    parser = argparse.ArgumentParser()
    parser.add_argument("parameters")
    return parser.parse_args()


def irange(start, stop, step=1):
    r = start
    while r <= stop:
        yield r
        r += step


def CONFIG_DESCRIPTIONS_TO_VARIABLES(parameters):
    DICTIONARY = {
        r"$\Lambda$":     parameters["target_growth_rate"],
        r"$d_{T}$":       parameters["death_rate"]["target"],
        r"$d_{S}$":       parameters["death_rate"]["sensitive"],
        r"$d_{R}$":       parameters["death_rate"]["resistant"],
        r"$d_{A}$":       parameters["death_rate"]["antibody"],
        r"$p_{S}$":       parameters["transition_mode_probability"]["sensitive"],
        r"$p_{R}$":       parameters["transition_mode_probability"]["resistant"],
        r"$\beta_{FS}$":  parameters["transition_rates"]["sensitive"]["free_viral"],
        r"$\beta_{CS}$":  parameters["transition_rates"]["sensitive"]["cell_to_cell"],
        r"$\beta_{FR}$":  parameters["transition_rates"]["resistant"]["free_viral"],
        r"$\beta_{CR}$":  parameters["transition_rates"]["resistant"]["cell_to_cell"],
        r"$\nu_{S}$":     parameters["antibody_effectivity"]["sensitive"],
        r"$\nu_{R}$":     parameters["antibody_effectivity"]["resistant"],
        r"$\mu$":         parameters["mutation_rate"],
        r"$N_{S}$":       parameters["burst_rate"]["sensitive"],
        r"$N_{R}$":       parameters["burst_rate"]["resistant"],
        r"$c_{S}$":       parameters["clearance_rate"]["sensitive"],
        r"$c_{R}$":       parameters["clearance_rate"]["resistant"],
        r"$q_{S}$":       parameters["antibody_absorption_rate"]["sensitive"],
        r"$q_{R}$":       parameters["antibody_absorption_rate"]["resistant"],
        r"$T_0$":         parameters["initial_values"]["target"],
        r"$I_{S0}$":      parameters["initial_values"]["sensitive_cells"],
        r"$I_{R0}$":      parameters["initial_values"]["resistant_cells"],
        r"$V_{S0}$":      parameters["initial_values"]["sensitive_virus"],
        r"$V_{R0}$":      parameters["initial_values"]["resistant_virus"],
        r"$A_0$":         parameters["initial_values"]["antibodies"],
    }

    return DICTIONARY


def plot_time_graphs(save_directory, step, system, text, args):
    plt.figure()

    plt.axes([0.20, 0.1, 0.75, 0.8], axisbg="white", frameon=True)

    limit = 1.1*max(list(system.T) + list(system.I_S) + list(system.I_R) + list(system.V_S) + list(system.V_R) + list(system.A))
    # limit = 1.1*max([max(system.T), max(system.I_S), max(system.I_R), max(system.V_S), max(system.V_R), max(system.A)])

    for index, text_line in enumerate(text):
        plt.text(-.25*system.t_f, limit*(1 - (0.05*index)), text_line)

    plt.ylim(-1., limit)
    plt.xlabel("Time")
    plt.ylabel("Densities")

    plt.plot(system.t, system.T, label=r"$T$", lw=2)
    plt.plot(system.t, system.I_S, label=r"$I_S$", lw=2)
    plt.plot(system.t, system.I_R, label=r"$I_R$", lw=2)
    plt.plot(system.t, system.V_S, label=r"$V_S$", lw=2)
    plt.plot(system.t, system.V_R, label=r"$V_R$", lw=2)
    plt.plot(system.t, system.A, label=r"$A$", lw=2)

    plt.legend(loc=0)

    filename = "%s/densities_%03d.png" % (save_directory, step)
    plt.savefig(filename, format="png", dpi=200)
    plt.close()
    garbage.collect()
    print GRAPH_SAVED % filename


def plot_graphs(save_directory, step, system, text, args):
    plot_time_graphs(save_directory, step, system, text, args)


class System():
    def __init__(self, parameters, t_f):
        self.t_f = t_f
        self.t = np.linspace(0, self.t_f, 100000)

        self.Lambda = parameters[r"$\Lambda$"]
        self.d_T = parameters[r"$d_{T}$"]
        self.d_S = parameters[r"$d_{S}$"]
        self.d_R = parameters[r"$d_{R}$"]
        self.d_A = parameters[r"$d_{A}$"]
        self.p_S = parameters[r"$p_{S}$"]
        self.p_R = parameters[r"$p_{R}$"]
        self.beta_FS = parameters[r"$\beta_{FS}$"]
        self.beta_CS = parameters[r"$\beta_{CS}$"]
        self.beta_FR = parameters[r"$\beta_{FR}$"]
        self.beta_CR = parameters[r"$\beta_{CR}$"]
        self.nu_S = parameters[r"$\nu_{S}$"]
        self.nu_R = parameters[r"$\nu_{R}$"]
        self.mu = parameters[r"$\mu$"]
        self.N_S = parameters[r"$N_{S}$"]
        self.N_R = parameters[r"$N_{R}$"]
        self.c_S = parameters[r"$c_{S}$"]
        self.c_R = parameters[r"$c_{R}$"]
        self.q_S = parameters[r"$q_{S}$"]
        self.q_R = parameters[r"$q_{R}$"]
        self.initial_values = [
            parameters[r"$T_0$"],
            parameters[r"$I_{S0}$"],
            parameters[r"$I_{R0}$"],
            parameters[r"$V_{S0}$"],
            parameters[r"$V_{R0}$"],
            parameters[r"$A_0$"]
        ]

        self.soln = odeint(self.f, self.initial_values, self.t)

        self.T = self.soln[:, 0]
        self.I_S = self.soln[:, 1]
        self.I_R = self.soln[:, 2]
        self.V_S = self.soln[:, 3]
        self.V_R = self.soln[:, 4]
        self.A = self.soln[:, 5]

    def f(self, y, t):
        T = y[0]
        I_S = y[1]
        I_R = y[2]
        V_S = y[3]
        V_R = y[4]
        A = y[5]

        next = [0]*6

        next[0] = self.Lambda - T*(self.d_T + self.p_S*(self.beta_FS/(1 + self.nu_S*A))*V_S + self.p_R*(self.beta_FR/(1 + self.nu_R*A))*V_R + (1 - self.p_S)*self.beta_CS*I_S + (1 - self.p_R)*self.beta_CR*I_R)
        next[1] = -self.d_S*I_S + T*(self.p_S*(self.beta_FS/(1 + self.nu_S*A))*V_S*(1 - self.mu) + (1 - self.p_S)*self.beta_CS*I_S*(1 - self.mu))
        next[2] = -self.d_R*I_R + T*(self.p_S*(self.beta_FS/(1 + self.nu_S*A))*V_S*self.mu + (1 - self.p_S)*self.beta_CS*I_S*self.mu + self.p_R*(self.beta_FR/(1 + self.nu_R*A))*V_R + (1 - self.p_R)*self.beta_CR*I_R)
        next[3] = self.p_S*self.N_S*self.d_S*I_S - self.c_S*V_S - self.p_S*self.beta_FS*V_S*T - self.q_S*V_S*A
        next[4] = self.p_R*self.N_R*self.d_R*I_R - self.c_R*V_R - self.p_R*self.beta_FR*V_R*T - self.q_R*V_R*A
        next[5] = A*(-self.d_A - self.q_S*V_S - self.q_R*V_R)

        return next


def main():
    # Gather the command line arguments
    args = PARSE_ARGS()

    # Load the parameters file (must be in json format)
    parameters = json.loads(open(args.parameters).read())

    # Get the current time to uniquely name a directory
    now = datetime.now()
    date_time_stamp = now.strftime(DATE_TIME_DIRECTORY_FORMAT)

    # Make the directory
    save_directory = "%s/../graphs/%s" % (DIRECTORY, date_time_stamp)
    os.system("mkdir -p %s" % save_directory)

    # Save the file that produced the graphs into the directory
    relevant_data_file = "%s/relevant_data.json" % save_directory
    pretty_data = json.dumps(parameters, indent=4, sort_keys=True)
    with open(relevant_data_file, "w") as relevant_data_file_object:
        relevant_data_file_object.write(pretty_data)

    # Refer to parameters by their variable names rather than their descriptions
    config_descriptions_to_variables = CONFIG_DESCRIPTIONS_TO_VARIABLES(parameters)

    # Read the final time from the config file
    t_f = parameters["final_time"]

    # Read the numbers of steps from the configuration file
    steps = int(parameters["steps"])

    # Determine the starting value, and the magnitude of the parameter shift
    starts_and_steps = {}
    for variable, dict_ in config_descriptions_to_variables.iteritems():
        starts_and_steps[variable] = {
            "start": dict_[0],
            "step": (dict_[1] - dict_[0])/steps
        }

    # Begin the loop
    for step in irange(0, steps):
        # Set current parameters based on the parameter shift
        current_parameters = {}
        for variable, dict_ in starts_and_steps.iteritems():
            current_parameters[variable] = dict_["start"] + (step*dict_["step"])

        # Perform the numerical integration
        system = System(current_parameters, t_f)

        # Get the text of all the changing parameters
        text = []
        for variable in current_parameters:
            if starts_and_steps[variable]["step"] != 0:
                text.append(variable + r"$=%.3f$" % current_parameters[variable])

        # Plot and save the graphs
        plot_graphs(save_directory, step, system, text, args)

    # For reference, print the date-time stamp
    print date_time_stamp
    print "\a\a\a"

if __name__ == "__main__":
    main()
