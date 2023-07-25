class Branch(BasePrimitive):
    """ Push events to HPQ  """
    
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = contest.pipeline_logger
        
    def _perform(self):
        self.action.args.extension = np.arange(1,9)
        
        for ext in self.action.args.extension:
            self.push_event('get_header', self.action.args.hdul[ext])
        return self.action.args