#
# User: dimkatsi91
# Date: December 22, 2020 | Covid Year 
# Simple example on how to plot a Pie Chart using matplotlib
#
# Update #1: Add bar chart example default & reverse orientation
#
import sys
from matplotlib import pyplot as plt

import numpy as np
# import pandas as pd


class Charts:
    """
        A class to plot a Pie Chart via matplotlib
        for a few Football teams & a metrics unit, like for
        example its people for the respective values
    """

    def __init__(self, teams, metrics, colors):
        """
            Constructor
        """
        self.teams = teams
        self.metrics = metrics
        self.colors = colors
        self.pre_check()

    def pre_check(self):
        """
            Check that all list have the same number of elements
            Mandatory condition | Other wise exit
        """
        if ( len(self.teams) is not len(self.metrics) ) or \
           ( len(self.teams) is not len(self.colors) ) or \
           ( len(self.metrics) is not len(self.colors) ):
           print("Pre-check condition not met.")
           print("Please check teams, metrics or colors lists for appropriate length and try again! Aborting ...")
           sys.exit(1)
        else:
            print("Pre-check condition valid. Continue procedure ...\n")

    def show_info(self):
        """
            Informational trivial function
        """
        print("\nPrinting-out the teams user has chosen ..")
        counter = 1
        for team in self.teams:
            print("Team %d : %s " %(counter, team))
            counter += 1
        #
        # Colors
        #
        print("\nPrinting-out the colors user has chosen ..")
        counter = 1
        for color in self.colors:
             print("Color %d : %s " %(counter, color))
             counter += 1
        #
        # Metrics for the respective teams
        #
        print("\nPrinting-out metrics of each flavor ..")
        counter = 0
        for metric in self.metrics:
            print("Team: %s | People: %d %%" %(self.teams[counter], metric))
            counter += 1

    def plot_pie_chart(self):
        """
            Plot the Pie Chart
        """
        plt.figure(figsize =(10, 7))
        explode = (0, 0.1, 0, 0, 0)
        plt.pie(self.metrics, explode=explode, labels=self.teams, colors=self.colors) # , startangle=90)
        plt.title("Greek Football Teams Sample Chart")
        plt.legend(self.teams, loc="best")
        plt.show()

    def list_to_tuple(self, input_list):
        """
            Input: a list
            Output: a tuple
        """
        return tuple(input_list)

    def bar_chart(self, align_option="default"):
        """
            Create a bar chart using Python's libraries (Numpy)
        """
        teams = self.list_to_tuple(self.teams)
        y_pos = np.arange(len(teams))
        metrics = self.metrics

        if align_option=="default":
            plt.bar(y_pos, metrics, align='center', alpha=0.5)
            plt.xticks(y_pos, teams)
            plt.ylabel('Percentage (%) --->')
            plt.xlabel('Teams ---> ')

        elif align_option=="reverse":
            plt.barh(y_pos, metrics, align='center', alpha=0.5)
            plt.yticks(y_pos, teams)
            plt.xlabel('Percentage (%) --->')
            plt.ylabel('Teams --->')

        else:
            print("Please choose valid option: [default | reverse] and try again!")
            sys.exit(1)
        
        plt.title('Greek Football Teams')
        plt.show()
    

def main():
    """
        main def()
    """
    #
    # Set a few flavors and the respective Kg for each of it
    # Last, set the colors to be used
    #
    my_teams = ["PAOK", "Olympiacos", "AEK", "Volos NFS", "PAO"]
    my_metrics = [15, 30, 20, 15, 20]
    my_colors = ["black", "red", "yellow", "navy", "green"]
    #
    # Create the object and call the function to plot the Pie Chart
    #
    my_chart = Charts(my_teams, my_metrics, my_colors)
    my_chart.show_info()
    my_chart.plot_pie_chart()
    my_chart.bar_chart()
    my_chart.bar_chart("reverse")

if __name__ == "__main__":
    main()
