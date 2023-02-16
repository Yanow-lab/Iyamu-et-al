__name__ = "ArrayTools"
__author__ = "Daniel Ferrer-Vinals"
__all__ = ["ArrayPlotter", "plot_alignment_array", "get_cmap", "get_colorbar" ]


import numpy as np
from biotite.sequence.graphics import LetterPlotter
from biotite.visualize import colors
from biotite.sequence.graphics.colorschemes import get_color_scheme

class ArrayPlotter(LetterPlotter):
    '''
    This SymbolPlotter maps data from epitope scannings on arrays 
    of overlapping peptides, onto proteins represented as sequencece aligment.
    Symbols are visualized as characters on a colored background box. The 
    color of a given box represent the recognition signal of the respective 
    20-mer peptide ending at that residue. The intensity of the color is 
    proportional to the fluorescence intensity recorded on the peptide array
    Fluorescence values are proportional to antibody recognition on 
    the peptide array.
    
    Parameters
    ----------
    axes : Axes
        A Matplotlib axes, that is used as plotting area.
    fl_score : numpy.ndarray 
        The ndarray to map fluorescence values to score residues.
        By default the normalized score is 1 for maximun recognition
        and 0 for non-recognition(no color).
    color_symbols : bool, optional
        If true, the symbols themselves are colored.
        If false, the symbols are black, and the boxes behind the
        symbols are colored.
    font_size : float, optional
        Font size of the sequence symbols.
    font_param : dict, optional
        Additional parameters that is given to the
        matplotlib.Text instance of each symbol.
    
    '''
    def __init__(self, axes, fl_score, color_symbols=False,
                 font_size=None, font_param=None):

        super().__init__(axes, color_symbols, font_size, font_param)
        if fl_score is not None:
            self.fl_score = fl_score #later develop signalarray.fl_score()
        else:
            self.fl_score = None        
        # Default colormap
        #self._cmap = _cmap
        self._cmap = self._generate_colormap(colors["dimorange"],
                                             self._color_symbols) # default= color_symbols: False
        
     
                             
    def get_color(self, alignment, column_i, seq_i):
        '''
        Get the color of a symbol at a specified position in the
        alignment.
        The symbol is specified as position in the alignment's trace
        (trace[pos_i, seq_i]).

        Parameters
        ----------
        alignment : Alignment
            The respective alignment.
        column_i : int
            The position index in the trace.
        seq_i : int
            The sequence index in the trace.
        Returns
        -------
        color : object
            A *Matplotlib* compatible color used for the background
            or at the specifed position
            '''
        index1 = alignment.trace[column_i, seq_i]
        if index1 == -1:
            spot_signal = 0
        else:
            spot_signal = self._get_signal(self.fl_score, column_i, seq_i)
           # spot_signal = 0.2  # from 0 to 0.9 it gives shades of orange
        return self._cmap(spot_signal)


    def _get_signal(self, fl_score, column_i, seq_i):
        if fl_score is None:
            signal = 0.01
        else:
            signal = fl_score[column_i, seq_i]
        return signal


            
    def plot_symbol(self, bbox, alignment, column_i, seq_i):
        '''
        This method is optimized to draw the symbols of an alignment while
        mapping the epitope scanning data throughout the entire sequence 
        alignment.
        '''
        from matplotlib.patches import Rectangle

        trace = alignment.trace
#        key3 = seq_i

        if trace[column_i, seq_i] != -1:
            key1 = alignment.sequences[1][trace[column_i, 1]]
            key2 = alignment.sequences[0][trace[column_i, 0]]
            if key1==key2:
#                key3 = seq_i
#                if key3 == 1:
                if seq_i == 1:
                    symbol = "*"
                else:
                    symbol = alignment.sequences[seq_i][trace[column_i, seq_i]]
            else:
                symbol = alignment.sequences[seq_i][trace[column_i, seq_i]]
        else:
            symbol = "-"
        color = self.get_color(alignment, column_i, seq_i) # where the color is comming from

        box = Rectangle(bbox.p0, bbox.width, bbox.height)
        self.axes.add_patch(box)
        text = self.axes.text(
            bbox.x0 + bbox.width/2, bbox.y0 + bbox.height/2,
            symbol, color="black", ha="center", va="center",
            size=self._font_size, **self._font_param)
        text.set_clip_on(True)

        if self._color_symbols:
            box.set_color("None")
            text.set_color(color)
        else:
            box.set_color(color) # here is where you color the box
        
    @staticmethod
    def _generate_colormap(color, to_black): # to_black is an alias for color_symbols
        from matplotlib.colors import ListedColormap, to_rgb

        color = to_rgb(color)
        if to_black: #this is always false as color_symbols = False default
            # From color to black
            cmap_val = np.stack(
                [np.interp(np.linspace(0, 1, 100), [0, 1], [color[i], 0])
                     for i in range(len(color))]
            ).transpose()
        else:
            # From white to color
            cmap_val = np.stack(
                [np.interp(np.linspace(0, 1, 100), [0, 1], [1, color[i]])
                    for i in range(len(color))]
            ).transpose()
        return ListedColormap(cmap_val)        




def plot_alignment_array(axes, alignment, symbols_per_line=50,
                                    show_numbers=False, number_size=None,
                                    number_functions=None,
                                    labels=None, label_size=None,
                                    show_line_position=False,
                                    spacing=1,
                                    color=None, cmap=None, fl_score= None,
                                    color_symbols=False, symbol_spacing=None,
                                    symbol_size=None, symbol_param=None):

    '''
    Plots a pairwise sequence alignment using an ArrayPloter instance.
    Higlights sequence recognition regions at the positions of the respective 
    score residue per alignment column.
    
    Parameters
    ----------
    axes : Axes
        A Matplotlib axes, that is used as plotting area.
    alignment : Alignment
        The pairwise sequence alignment to be plotted.
    symbol_plotter : SymbolPlotter
        Instance of ArrayPlotter. Defines how the symbols are drawn
        in the alignment.
    symbols_per_line : int, optional
        The amount of alignment columns that are diplayed per line.
    show_numbers : bool, optional
        If true, the sequence position of the symbols in the last
        alignment column of a line is shown on the right side of the
        plot.
        If the last symbol is a gap, the position of the last actual
        symbol before this gap is taken.
        If the first symbol did not occur up to this point,
        no number is shown for this line.
        By default the first symbol of a sequence has the position 1,
        but this behavior can be changed using the `number_functions`
        parameter.
    number_size : float, optional
        The font size of the position numbers
    number_functions : list of [(None or Callable(int -> int)], optional
        By default the position of the first symbol in a sequence is 1,
        i.e. the sequence position is the sequence index incremented by
        1.
        The behavior can be changed with this parameter:
        If supplied, the length of the list must match the number of
        sequences in the alignment.
        Every entry is a function that maps a sequence index (*int*) to
        a sequence position (*int*) for the respective sequence.
        A `None` entry means, that the default numbering is applied
        for the sequence.
    labels : list of str, optional
        The sequence labels.
        Must be the same size and order as the sequences in the
        alignment.
    label_size : float, optional
        Font size of the labels
    show_line_position : bool, optional
        If true the position within a line is plotted below the
        alignment.
    spacing : float, optional
        The spacing between the alignment lines. 1.0 means that the size
        is equal to the size of a symbol box.
    color : tuple or str, optional
        A *Matplotlib* compatible color.
    cmap : Colormap or str, optional
        The boxes (or symbols, if `color_symbols` is set) are
        colored based on the normalized intensity value on the
        given *Matplotlib* Colormap.
    fl_score : numpy.ndarray 
        The ndarray to map fluorescence values to score residues.
        By default the normalized score is 1 for maximun recognition
        and 0 for non-recognition(no color).
    color_symbols : bool, optional
        If true, the symbols themselves are colored.
        If false, the symbols are black, and the boxes behind the
        symbols are colored.
    symbol_size : float, optional
        Font size of the sequence symbols.
    symbol_param : dict
        Additional parameters that is given to the
        :class:`matplotlib.Text` instance of each symbol.
    symbol_spacing : int, optional
        А space is placed between each number of elements desired
        by variable.
    Notes
    -----
    A '*' represents a sequence match on the alignment 
    A '-' represents a sequence gap on the alignment 

    '''
    symbol_plotter = ArrayPlotter(
        axes, fl_score = fl_score, font_size = symbol_size, font_param = symbol_param,
        color_symbols = color_symbols
    )
    
    if color is not None or cmap is not None:
        symbol_plotter.set_color(color=color, cmap=cmap)
        
    from biotite.sequence.graphics import plot_alignment
    plot_alignment(
        axes=axes, alignment=alignment, symbol_plotter=symbol_plotter,
        symbols_per_line=symbols_per_line,
        show_numbers=show_numbers, number_size=number_size,
        number_functions=number_functions,
        labels=labels, label_size=label_size,
        show_line_position=show_line_position,
        spacing=spacing, symbol_spacing=symbol_spacing
    )  
    
    
def get_cmap(axes, signal_map):
    '''
    Access the Colormap atribute of an ArrayPlotter instantiated: 
    ArrayPlotter(axes, signal_map). 
    ----------
    axes: Axes
        A Matplotlib axes to plot the Colormap
    signal_map: numpy.ndarray
        A numpy.ndarray that maps peptide signal score to the corresponding 
        position of the score residue on the sequence alignment
  
    Returns
    ----------
    A matplotlib.colors.Colormap object, used to convert data from the signal_map
    to the RGBA colors 
    '''
    ax = axes
    smap = signal_map
    plotter = ArrayPlotter(ax, smap)
    return plotter._cmap


def get_colorbar(axes, array1, array2, colormap, transform ='linear',  orient =None,
                title=None):
    '''
    Generate the a Colorbar object according to the specified scaling of the data
    ----------
    axes: Axes
        A Matplotlib axes were the Colorbar will reside
    array1: pandas.DataFrame
        A dataframe that contains the peptide array scan data(antigen-1) with the desired
        trasnfomation method
    array2: pandas.DataFrame
        A dataframe that contains the peptide array scan data(antigen-2) with the desired
        trasnfomation method
    transform: str. Optional
        One of: 'linear', 'sqrt', 'cubic', 'log'(Default: 'linear').Instantiate subclasses
        of matplotlib.cm.ScalarMappable to scale the colors on the desired transformation
    orient: 'vertical' or 'horizontal'
        Orientation of the Colorbar
    title: str
        The label on the Colorbar's long axis
    
    Returns
    ----------
    A matplotlib.Colorbar instance. Draw a colorbar in an existing axes
    '''
    import matplotlib as mpl
    
    df1 = array1
    df2 = array2
    cmp = colormap
    method = transform
    ax = axes
    orientation = orient
    label = title
    
    # custom Formtatter for tick labels on the colorbar
    def fmt(x, pos):
        a, b = '{:.1e}'.format(x).split('e')
        b = int(b)
        return r'${}\cdot10^{{{}}}$'.format(a, b)
    
    if method == 'linear':
        vmiA = df1['combined_signal'].min()
        vmiB = df2['combined_signal'].min()
        vmxA = df1['combined_signal'].max()
        vmxB = df2['combined_signal'].max()
        # Colormap normalization:
        norm = mpl.colors.PowerNorm(gamma = 1.0, 
                                    vmin = min(vmiA,vmiB), vmax = max(vmxA,vmxB))

    elif method == 'sqrt':
        vmiA = df1['combined_signal'].min()
        vmiB = df2['combined_signal'].min()
        vmxA = df1['combined_signal'].max()
        vmxB = df2['combined_signal'].max()
        # Colormap normalization:
        norm = mpl.colors.PowerNorm(gamma = 0.5, 
                                    vmin = min(vmiA,vmiB), vmax = max(vmxA,vmxB))
    elif method == 'cubic':
        vmiA = df1['combined_signal'].min()
        vmiB = df2['combined_signal'].min()
        vmxA = df1['combined_signal'].max()
        vmxB = df2['combined_signal'].max()
        # Colormap normalization:
        norm = mpl.colors.PowerNorm(gamma = 0.33, 
                                    vmin = min(vmiA,vmiB), vmax = max(vmxA,vmxB))
    elif method == 'log':
        vmiA = df1['combined_signal'].min()
        vmiB = df2['combined_signal'].min()
        vmxA = df1['combined_signal'].max()
        vmxB = df2['combined_signal'].max()
        # Colormap normalization:
        norm = mpl.colors.LogNorm(vmin = 1 if vmiA==0 or vmiB==0 else min(vmiA,vmiB),
                                  vmax = max(vmxA,vmxB))
        
    fig = mpl.pyplot.figure()        
    return fig.colorbar(mpl.cm.ScalarMappable(norm = norm, cmap = cmp),
                 cax = ax, orientation = orientation, label = label,
                 format = mpl.ticker.FuncFormatter(fmt))
    
    