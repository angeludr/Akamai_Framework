class RemoveOverscan(BasePrimitive):
    """ Remove overscan regions from CCD images """
    
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
        
    def _perform(self):
        self.action.args.data = self.action.args.data[:, self.action.args.precol: self.action.args.width - self.action.args.postpix]
        self.action.args.handled += 1
        
        self.action.args.alldata = []
        self.action.args.alldata.append(self.action.args.data)
        return self.action.args