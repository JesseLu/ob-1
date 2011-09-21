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
        """ Plot data. Input parameter data should be a numpy array of either
            one or two dimensions.

        """

        if data.ndim == 1:
            self.g.plot(data)
        elif data.ndim == 2:
            self.g.plot(Gnuplot.GridData(data, range(data.shape[0]), 
                    range(data.shape[1]), binary=0, with_='image'))
        else:
            raise Exception("Data must be a 1 or 2-dimensional array,")

    def save(self, filename):
        self.g.hardcopy(filename, enhanced=1, color=1)
        


def create(debug_option=0):
    """ Create a gnuplotter class instance.

    """
    return Plotter(debug_option=debug_option)


