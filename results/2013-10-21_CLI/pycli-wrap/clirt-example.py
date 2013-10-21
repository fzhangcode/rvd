#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

#check operating system
if os.name == 'posix':
    from lib import CLIposix as CLIwrap
elif os.name == 'nt':
    from lib import CLInt as CLIwrap
else:
    print "Error: current OS not supported ;)"
    sys.exit(0)
      
cli = CLIwrap.CLImain()

configure_dict = {
    'PromptStr': "myCLI> ",
    'InitialTasks': cli.InitialTasks,
    'PressEnterTasks': cli.PressEnterTasks,
    'CleanUpTasks': cli.CleanUpTasks,
    'InfiniteLoopTasks': cli.InfiniteLoopTasks,
    'ProcessCLI': cli.ProcessCLI,
    'EnterKeys': cli.EnterKeys,
    'ExitWords': cli.ExitWords,
    'ExitKeys': cli.ExitKeys,
}

cli.configure(**(configure_dict))
cli.run()