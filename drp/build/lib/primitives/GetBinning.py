class GetBinning(BasePrimitive):
    """ Get binning dimensions """
    
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
    
    def _perform(self):
        self.action.args.binning = self.action.args.hdr0['BINNING'].split(',')
        return self.section.args