class OpenFits(BasePrimitive):
    """ Opens fits file """
    
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
        
    def _perform(self):
        self.action.args.hdul = fits.open(self.action.args.filename)
        self.action.args.handled = 0
        return self.action.args