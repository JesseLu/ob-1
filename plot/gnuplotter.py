import Gnuplot, Gnuplot.funcutils


class Plotter:
    """ A class to manage plotting.

    """

    def __init__(self, debug_option=0):
        """ Create Gnuplot instance (using the Gnuplot for python package).
            
            Use debug_option=1 to enable debugging output, default is False.

        """
        self.g = Gnuplot.Gnuplot(debug=debug_option)
        self.g('set data style linespoints') # Sets option for 1D plotting.


    def plot(self, data):
        self.g.plot(data)

        


def create(debug_option=0):
    """ Create a gnuplotter class instance.

    """
    return Plotter(debug_option=debug_option)


