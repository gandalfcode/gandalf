class Data:
    '''Generic class to communicate the data with the plotting part.
    Simple container with x_data, y_data and render_data members,
    which are supposed to be numpy arrays.
    '''
    def __init__(self, x_data, y_data, render_data = None):
        self.x_data = x_data
        self.y_data = y_data
        self.render_data = render_data