class PlotBias(BasePrimitive):
    """ Plot bias' of each CCD """
    
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
        
    def _perform(self):
        plt.plot(self.action.args.bias)
        plt.xlim(xmin=-500, xmax=4500)
        plt.ylim(ymin=800, ymax=1800)
        plt.xlabel('PIXEL')
        plt.ylabel('#COUNTS')
        plt.show()