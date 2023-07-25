class GetOverscan(BasePrimitive):
    """ Get the parameters of the overscan region of each detector """
    
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
    
    def _perform(self):
        self.action.args.precol   = int(self.action.args.hdr0['PRECOL'])   // int(self.action.args.binning[0])
        self.action.args.postpix  = int(self.action.args.hdr0['PPOSTPIX']) // int(self.action.args.binning[0])
        self.action.args.preline  = int(self.action.args.hdr0['PRELINE'])  // int(self.action.args.binning[1])
        self.action.args.postline = int(self.action.args.hdr0['POSTLINE']) // int(self.action.args.binning[1])
        return self.action.args