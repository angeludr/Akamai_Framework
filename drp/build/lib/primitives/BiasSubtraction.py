class BiasSubtraction(BasePrimitive):
    """ Subtract the bias from the detector data """
    
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
        
    def _perform(self):
        self.action.args.data = self.action.args.data - self.action.args.bias[:, None]
        return self.action.args