class DetectorMosaic(BasePrimitive):
    """ Creating mosaic of DEIMOS 8 CCDs """
    
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
        
    def _perform(self):
        self.action.args.norm = ImageNormalize(self.action.args.data, 
                              self.action.args.vmin == self.action.args.vmin, 
                              self.action.args.vmax == self.action.args.vmax)
        fig = plt.figure(frameon=False)
        ax = fig.add_axes([0,0,1,1])
        plt.axis('off')
        plt.imshow(self.action.args.data, cmap='gray', norm=self.action.args.norm)
        plt.savefig(self.action.args.filename, + '.jpg', dpi=300)