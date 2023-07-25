class GetHeader(BasePrimitive):
    """ Get header of fits file """
    
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
        
    def _perform(self):
        hdr0 = self.action.args.hdul[self.action.args.extension].header
        return self.action.args