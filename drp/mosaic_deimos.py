from keckdrpframework.core.framework import Framework
from keckdrpframework.config.framework_config import ConfigClass
from keckdrpframework.models.arguments import Arguments
from keckdrpframework.utils.drpf_logger import getLogger
from primitives import *

from pipeline import Akamai_Pipeline

import argparse
import sys
import traceback



def _parse_arguments(in_args: list) -> argparse.Namespace:
    description = "DEIMOS Mosaicing pipeline CLI"

    # this is a simple case where we provide a frame and a configuration file
    parser = argparse.ArgumentParser(prog=f"{in_args[0]}",
                                     description=description)

    parser.add_argument('-f', '--frames', nargs='*', type=str,
                        help='input image files (full path, list ok)',
                        default=None)


    # after ingesting the files,
    # do we want to continue monitoring the directory?

    out_args = parser.parse_args(in_args[1:])
    return out_args

def main():

    # def process_subset(in_subset):
    #     for in_frame in in_subset.index:
    #         arguments = Arguments(name=in_frame)
    #         framework.append_event('next_file', arguments, recurrent=True)
    
    # def process_list(in_list):
    #     for in_frame in in_list:
    #         arguments = Arguments(name=in_frame)
    #         framework.append_event('next_file', arguments, recurrent=True)

    args = _parse_arguments(sys.argv)

    framework_config_file = "./configs/framework.cfg"


    try:
        framework = Framework(Akamai_Pipeline, framework_config_file)
    except Exception as e:
        print("Failed to initialize framework, exiting ...", e)
        traceback.print_exc()
        sys.exit(1)
    framework.context.pipeline_logger = getLogger(name="Akamai_Logger")
    framework.logger = getLogger('../configs/logger.ini', name="DRPF")

    # if args.infiles is not None:
    #     framework.config.file_type = args.infiles

    # # start the bokeh server is requested by the configuration parameters
    # if framework.config.instrument.enable_bokeh is True:
    #     if check_running_process(process='bokeh') is False:
    #         subprocess.Popen('bokeh serve', shell=True)
    #         # --session-ids=unsigned --session-token-expiration=86400',
    #         # shell=True)
    #         time.sleep(5)
    #     # subprocess.Popen('open http://localhost:5006?bokeh-session-id=kcwi',
    #     # shell=True)

    framework.logger.info("Framework initialized")

    # add a start_bokeh event to the processing queue,
    # # if requested by the configuration parameters
    # if framework.config.instrument.enable_bokeh:
    #     framework.append_event('start_bokeh', None)

    # single frame processing
    if args.frames:
        framework.append_event('open_fits_file', None)
        # framework.open_fits_file(None, args.frames, False)

    framework.start(False, False, False, False)


if __name__ == "__main__":
    main()