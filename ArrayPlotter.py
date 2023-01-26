#__name__ = "array..."
__author__ = "Daniel Ferrer-Vinals"
__all__ = ["ArrayPlotter", "plot_alignment_array" ]


import numpy as np
from biotite.sequence.graphics import LetterPlotter
from biotite.visualize import colors
from biotite.sequence.graphics.colorschemes import get_color_scheme

class ArrayPlotter(LetterPlotter):
    """
    This :class:`SymbolPlotter` plots peptide recognition data  
    for two protein sequences represented as global aligment 
    """
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
        index1 = alignment.trace[column_i, seq_i]
        if index1 == -1:
            spot_signal = 0
        else:
            spot_signal = self._get_signal(self.fl_score, column_i, seq_i)
#            spot_signal = 0.2                                          # from 0 to 0.9 it gives shades of green
        return self._cmap(spot_signal)


    def _get_signal(self, fl_score, column_i, seq_i):
        if fl_score is None:
            signal = 0.01
        else:
            signal = fl_score[column_i, seq_i]
        return signal


            
    def plot_symbol(self, bbox, alignment, column_i, seq_i):
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
                [np.interp(np.linspace(0, 1, 10), [0, 1], [1, color[i]])
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

    """
    Plot a pairwise sequence alignment highlighting
    the sequence recognition regions per alignment column
    """
    symbol_plotter = ArrayPlotter(
        axes, fl_score= fl_score, font_size=symbol_size, font_param=symbol_param,
        color_symbols=color_symbols
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