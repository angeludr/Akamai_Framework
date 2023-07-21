class StackRows(BasePrimitive):
    """ Stacking top and bottom rows of mosaic together and rotating the final image """
    
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
        
    def _perform(self):
        self.action.args.fulldata = np.concatenate((self.action.args.top, self.action.args.bottom), axis=0)
        self.action.args.fulldata = np.rot90(self.action.args.fulldata)
        return self.action.args